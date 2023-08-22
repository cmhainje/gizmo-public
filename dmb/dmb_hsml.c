#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! \file dmb_hsml.c
 *  \brief kernel length determination for dark matter-baryon interactions
 *
 *  This file contains a loop modeled on the gas density computation which
 *    determines softening lengths (and appropriate correction terms)
 *    for the nearest DM particles about a gas cell and vice versa. It then
 *    computes the momentum and heat transfer rates between DM and gas.
 */
/*
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#ifdef DM_DMB

#define CORE_FUNCTION_NAME dmb_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME dmb_particle2in    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME dmb_out2particle  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(dmb_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

double DMB_MAX_ALLOWED_GAS_HSML = 10.0;

int dmb_BITFLAG(short int particle_type)
{
    // if (particle_type == 0) {
    //     return 2; /* 2^1: gas sees DM */
    // }
    if (particle_type == 1) {
        return 1; /* 2^0: DM sees gas */
    }
    return 0;
}

/* this structure defines the variables that need to be sent -from- the 'searching' element */
static struct INPUT_STRUCT_NAME
{
    MyDouble Pos[3];
    MyDouble Vel[3];
    MyDouble Mass;
    MyDouble Temp;
    // MyFloat Hsml;
    MyFloat AGS_Hsml;
    MyIDType ID;
    int NodeList[NODELISTLENGTH];
    int Type;

}
*DATAIN_NAME, *DATAGET_NAME;

/* this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
void dmb_particle2in(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k; for (k = 0; k < 3; k++) {
        in->Pos[k] = P[i].Pos[k];
        in->Vel[k] = P[i].Vel[k];
    }
    in->Mass = P[i].Mass;
    in->Temp = P[i].DMB_MyTemp;
    // in->Hsml = P[i].DMB_Hsml;
    in->AGS_Hsml = (P[i].Type == 0 ? PPP[i].Hsml : PPP[i].AGS_Hsml);
    in->ID = P[i].ID;
    in->Type = P[i].Type;
}


/* this structure defines the variables that need to be sent -back to- the 'searching' element */
static struct OUTPUT_STRUCT_NAME
{
    MyLongDouble Ngb;
    int NgbInt;
    MyLongDouble V[3];
    MyLongDouble Density;
    MyLongDouble Temp;
    MyLongDouble GasMass;
}
*DATARESULT_NAME, *DATAOUT_NAME;

/* this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
void dmb_out2particle(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    ASSIGN_ADD(P[i].DMB_NumNgb, out->Ngb, mode);
    ASSIGN_ADD(P[i].DMB_NgbInt, out->NgbInt, mode);
    int k; for (k = 0; k < 3; k++) { ASSIGN_ADD(P[i].DMB_V[k], out->V[k], mode); }
    ASSIGN_ADD(P[i].DMB_Density, out->Density, mode);
    ASSIGN_ADD(P[i].DMB_Temperature, out->Temp, mode);
    ASSIGN_ADD(P[i].DMB_GasMass, out->GasMass, mode);
}


/* routine to determine if we need to use dmb_evaluate to calculate Hsml */
int dmb_isactive(int i);
int dmb_isactive(int i)
{
    if(P[i].TimeBin < 0) return 0;
    if(P[i].Type != 1) return 0; /* only DM particles */
    if(P[i].Mass <= 0) return 0;
    return 1;
}

int dmb_isactive_gas(int i);
int dmb_isactive_gas(int i)
{
    if(P[i].TimeBin < 0) return 0;
    if(P[i].Type != 0) return 0; /* only gas particles */
    if(P[i].Mass <= 0) return 0;
    return 1;
}

struct kernel_pair
{
    double dp[3], dv[3], r, wk_i, wk_j, dwk_i, dwk_j, h_i, hinv_i, hinv3_i, hinv4_i, h_j, hinv_j, hinv3_j, hinv4_j;
};

/*! This function represents the core of the density computation. The target particle may either be local, or reside in the communication buffer. */
/*!   -- this subroutine contains no writes to shared memory -- */
int dmb_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, n, k;
    struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    if(mode == 0)
        INPUTFUNCTION_NAME(&local, target, loop_iteration);
    else
        local = DATAGET_NAME[target];

    int bitflag = dmb_BITFLAG(local.Type);

    struct kernel_pair kernel;
    kernel.h_i = local.AGS_Hsml;
    kernel_hinv(kernel.h_i, &kernel.hinv_i, &kernel.hinv3_i, &kernel.hinv4_i);
    double r2, u_i, u_j, g_ij, h_ratio;
    double search_len = local.AGS_Hsml + DMBMaxGasHsml;

    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, search_len, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, bitflag);
            if(numngb_inbox < 0) {return -2;}
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if((P[j].Mass <= 0) || (PPP[j].Hsml <= 0)) continue;

                for (k = 0; k < 3; k++) { kernel.dp[k] = local.Pos[k] - P[j].Pos[k]; }
                NEAREST_XYZ(kernel.dp[0], kernel.dp[1], kernel.dp[2], 1);
                r2 = kernel.dp[0]*kernel.dp[0] + kernel.dp[1]*kernel.dp[1] + kernel.dp[2]*kernel.dp[2];
                if (r2 <= 0) continue;
                kernel.r = sqrt(r2);
                kernel.h_j = PPP[j].Hsml;
                kernel_hinv(kernel.h_j, &kernel.hinv_j, &kernel.hinv3_j, &kernel.hinv4_j);

                if (kernel.r > kernel.h_i + kernel.h_j) continue;


                // compute the overlap
                h_ratio = kernel.h_i / kernel.h_j;
                g_ij = 0.0;
                if (h_ratio >= 2) {
                    u_i = kernel.r * kernel.hinv_i;
                    if (u_i < 1) {
                        kernel_main(u_i, kernel.hinv3_i, kernel.hinv4_i, &kernel.wk_i, &kernel.dwk_i, 0);
                        g_ij = kernel.wk_i;
                    }
                }
                else if (h_ratio <= 0.5) {
                    u_j = kernel.r * kernel.hinv_j;
                    if (u_j < 1) {
                        kernel_main(u_j, kernel.hinv3_j, kernel.hinv4_j, &kernel.wk_j, &kernel.dwk_j, 0);
                        g_ij = kernel.wk_j;
                    }
                }
                else {
                    g_ij = g_geo(kernel.r * 2 / (kernel.h_i + kernel.h_j));
                }


                /* DMB core computations: density, velocity, temperature, numngb for both particles */
                if (g_ij != 0) {
                    // update DM particle
                    out.Density += g_ij * P[j].Mass;
                    for (k = 0; k < 3; k++) { out.V[k] += g_ij * P[j].Vel[k]; }
                    out.Temp += g_ij * P[j].DMB_MyTemp;
                    out.GasMass += g_ij * P[j].DMB_MyMass;
                    out.Ngb += g_ij;
                    out.NgbInt++;

                    // update gas particle; taking care with writes to shared memory
                    #pragma omp atomic update
                    P[j].DMB_Density += g_ij * local.Mass;

                    for (k = 0; k < 3; k++) {
                        #pragma omp atomic update
                        P[j].DMB_V[k] += g_ij * local.Vel[k];
                    }

                    #pragma omp atomic update
                    P[j].DMB_Temperature += g_ij * local.Temp;

                    #pragma omp atomic update
                    P[j].DMB_NumNgb += g_ij;

                    #pragma omp atomic update
                    P[j].DMB_NgbInt += 1;
                }
            } // numngb_inbox loop
        } // while(startnode)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    }
    if(mode == 0) {OUTPUTFUNCTION_NAME(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
}


void dmb_calc(void)
{
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second();
    PRINT_STATUS(" ..entering DMB interaction calculation\n");

    /* initialize anything we need to about the active particles before their loop */
    DMBMaxGasHsml = 0.0;
    int i; for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        if (dmb_isactive(i) | dmb_isactive_gas(i)) {
            P[i].DMB_NumNgb = 0;

            // compute temperatures once per particle
            if (P[i].Type == 0) {
                double rho = SphP[i].Density * All.cf_a3inv;
                double u = SphP[i].InternalEnergyPred;
                double mu=1, ne=1, nh0=0, nHe0, nHepp, nhp, nHeII;
                double T = ThermalProperties(u, rho, i, &mu, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp) * BOLTZMANN_CGS;
                P[i].DMB_MyTemp = T;
                P[i].DMB_MyMass = mu * PROTONMASS_CGS;

                DMBMaxGasHsml = DMAX(DMBMaxGasHsml, PPP[i].Hsml);
            }
            else if (P[i].Type == 1) {
                P[i].DMB_MyTemp = temperature_DM(P[i].AGS_VelDisp);
                if (isnan(P[i].DMB_MyTemp)) {
                    printf("DM temperature is nan\n");
                    printf("  AGS_VelDisp = %.3e\n", P[i].AGS_VelDisp);
                }
            }
        }
    }
    DMBMaxGasHsml = DMIN(DMBMaxGasHsml, DMB_MAX_ALLOWED_GAS_HSML);

    /* run the nearest neighbors loops */
    #include "../system/code_block_xchange_perform_ops_malloc.h"
    #include "../system/code_block_xchange_perform_ops.h"
    #include "../system/code_block_xchange_perform_ops_demalloc.h"

    /* final operations on the results */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(dmb_isactive(i) | dmb_isactive_gas(i))
        {
            if (P[i].DMB_NumNgb > 0) {
                int k; for (k = 0; k < 3; k++) { P[i].DMB_V[k] /= P[i].DMB_NumNgb; }
                P[i].DMB_Temperature /= P[i].DMB_NumNgb;
                if (P[i].Type == 1) { P[i].DMB_GasMass /= P[i].DMB_NumNgb; }
            } else {
                int k; for (k = 0; k < 3; k++) { P[i].DMB_V[k] = 0; }
                P[i].DMB_Density = 0;
                P[i].DMB_Temperature = 0;
                if (P[i].Type == 1) { P[i].DMB_GasMass = 0; }
            }

            // now all the ingredients are known -> compute and apply exchange rates!
            double acc[3], dUdt; compute_exch_rates(i, acc, &dUdt);
            int k; for (k = 0; k < 3; k++) { P[i].GravAccel[k] += acc[k]; }
            if (P[i].Type == 0) { SphP[i].DtInternalEnergy += dUdt; }
            else                { P[i].DMB_InternalEnergy += dUdt; }
        }
    }

    /* collect timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp; CPU_Step[CPU_AGSDENSWAIT] += timewait;
    CPU_Step[CPU_AGSDENSCOMM] += timecomm; CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */

#endif
