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

int dmb_BITFLAG(short int particle_type)
{
    if (particle_type == 0) {
        return 2; /* 2^1: gas sees DM */
    }
    else if (particle_type == 1) {
        return 1; /* 2^0: DM sees gas */
    }
    return 0;
}

/* this structure defines the variables that need to be sent -from- the 'searching' element */
static struct INPUT_STRUCT_NAME
{
    MyDouble Pos[3];
    MyFloat Hsml;
    int NodeList[NODELISTLENGTH];
    int Type;
}
*DATAIN_NAME, *DATAGET_NAME;

/* this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
void dmb_particle2in(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k; for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k];}
    in->Hsml = P[i].DMB_Hsml;
    in->Type = P[i].Type;
}


/* this structure defines the variables that need to be sent -back to- the 'searching' element */
static struct OUTPUT_STRUCT_NAME
{
    MyLongDouble Ngb;
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
    if(P[i].Type > 1) return 0; /* only gas, DM particles */
    if(P[i].Mass <= 0) return 0;
    return 1;
}


/*! This function represents the core of the density computation. The target particle may either be local, or reside in the communication buffer. */
/*!   -- this subroutine contains no writes to shared memory -- */
int dmb_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, n;
    struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    if(mode == 0)
        INPUTFUNCTION_NAME(&local, target, loop_iteration);
    else
        local = DATAGET_NAME[target];

    int bitflag = dmb_BITFLAG(local.Type);

    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, bitflag);
            if(numngb_inbox < 0) {return -2;}
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if(P[j].Mass <= 0) continue;
                int k; for (k = 0; k < 3; k++) { out.V[k] += P[j].Vel[k]; }
                if (P[j].Type == 0) {
                    out.Density += SphP[j].Density;
                    out.Temp += P[j].DMB_MyTemp;
                    out.GasMass += P[j].DMB_MyMass;
                }
                else if (P[j].Type == 1) {
                    out.Density += P[j].AGS_Density;
                    out.Temp += P[j].DMB_MyTemp;
                }
                out.Ngb++;
            } // numngb_inbox loop
        } // while(startnode)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    }
    if(mode == 0) {OUTPUTFUNCTION_NAME(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
}


void dmb_calc(void)
{
    /* initialize variables used below, in particlar the structures we need to call throughout the iteration */
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second();
    MyFloat *Left, *Right; double desnumngb=64, desnumngbdev=48; long long ntot; 
    int i, npleft, iter=0, redo_particle;
    Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
    Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));

    /* initialize anything we need to about the active particles before their loop */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        if (dmb_isactive(i)) {
            P[i].DMB_NumNgb = 0;
            Left[i] = Right[i] = 0;

            // compute temperatures once per active particle
            if (P[i].Type == 0) {
                double rho = SphP[i].Density * All.cf_a3inv;
                double u = SphP[i].InternalEnergyPred;
                double mu=1, ne=1, nh0=0, nHe0, nHepp, nhp, nHeII;
                double T = ThermalProperties(u, rho, i, &mu, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp) * BOLTZMANN_CGS;
                P[i].DMB_MyTemp = T;
                P[i].DMB_MyMass = mu * PROTONMASS_CGS;
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
    
    /* allocate buffers to arrange communication */
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
    do
    {
        #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */

        /* do check on whether we have enough neighbors, and iterate for density-hsml solution */
        double tstart = my_second(), tend;
        for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
        {
            if(dmb_isactive(i))
            {
                redo_particle = 0; /* now check whether we have enough neighbours, and are below the maximum search radius */
                double maxsoft = DMIN(All.MaxHsml, 10.0*PPP[i].AGS_Hsml);

                if(((P[i].DMB_NumNgb < desnumngb - desnumngbdev) || (P[i].DMB_NumNgb > (desnumngb + desnumngbdev)))
                   && (Right[i]-Left[i] > 0.001*Left[i] || Left[i]==0 || Right[i]==0))
                {
                    redo_particle = 1;
                }
                if(P[i].DMB_Hsml >= maxsoft)
                {
                    P[i].DMB_Hsml = maxsoft;
                    redo_particle = 0;
                }

                if(redo_particle)
                {
                    /* need to redo this particle */
                    npleft++;
                    
                    if(P[i].DMB_NumNgb < desnumngb-desnumngbdev) {Left[i]=DMAX(P[i].DMB_Hsml, Left[i]);}
                    if(P[i].DMB_NumNgb > desnumngb+desnumngbdev) {if(Right[i]>0) {Right[i]=DMIN(Right[i],P[i].DMB_Hsml);} else {Right[i]=P[i].DMB_Hsml;}}
                    
                    if(iter >= MAXITER - 10)
                    {
                        PRINT_WARNING("DMB: i=%d task=%d ID=%llu Type=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
                         i, ThisTask, (unsigned long long) P[i].ID, P[i].Type, P[i].DMB_Hsml, Left[i], Right[i], (float) P[i].DMB_NumNgb, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                    }
                    
                    // right/left define upper/lower bounds from previous iterations
                    if(Right[i] > 0 && Left[i] > 0)
                    {
                        // geometric interpolation between right/left //
                        if(P[i].DMB_NumNgb > 1)
                        {
                            P[i].DMB_Hsml *= pow( desnumngb / P[i].DMB_NumNgb , 1./NUMDIMS );
                        } else {
                            P[i].DMB_Hsml *= 2.0;
                        }
                        if((P[i].DMB_Hsml<Right[i])&&(P[i].DMB_Hsml>Left[i]))
                        {
                            P[i].DMB_Hsml = pow(P[i].DMB_Hsml*P[i].DMB_Hsml*P[i].DMB_Hsml*P[i].DMB_Hsml * Left[i]*Right[i] , 1./6.);
                        } else {
                            if(P[i].DMB_Hsml>Right[i]) P[i].DMB_Hsml=Right[i];
                            if(P[i].DMB_Hsml<Left[i]) P[i].DMB_Hsml=Left[i];
                            P[i].DMB_Hsml = pow(P[i].DMB_Hsml * Left[i] * Right[i] , 1.0/3.0);
                        }
                    }
                    else
                    {
                        if(Right[i] == 0 && Left[i] == 0)
                        {
                            char buf[1000]; sprintf(buf, "DMB: Right[i] == 0 && Left[i] == 0 && P[i].DMB_Hsml=%g\n", P[i].DMB_Hsml); terminate(buf);
                        }
                        double fac;
                        if(Right[i] == 0 && Left[i] > 0)
                        {
                            if(P[i].DMB_NumNgb > 1) {fac = log( desnumngb / P[i].DMB_NumNgb ) / NUMDIMS;} else {fac=1.4;}
                            if((P[i].DMB_NumNgb < 2*desnumngb)&&(P[i].DMB_NumNgb > 0.1*desnumngb)) {P[i].DMB_Hsml *= exp(fac);} else {P[i].DMB_Hsml *= 1.26;}
                        }
                        
                        if(Right[i] > 0 && Left[i] == 0)
                        {
                            if(P[i].DMB_NumNgb > 1) {fac = log( desnumngb / P[i].DMB_NumNgb ) / NUMDIMS;} else {fac=1.4;}
                            fac = DMAX(fac,-1.535);
                            if((P[i].DMB_NumNgb < 2*desnumngb)&&(P[i].DMB_NumNgb > 0.1*desnumngb)) {P[i].DMB_Hsml *= exp(fac);} else {P[i].DMB_Hsml /= 1.26;}
                        }
                    }
                }
                else {P[i].TimeBin = -P[i].TimeBin - 1;}	/* Mark as inactive */
            } //  if(dmb_isactive(i))
        } // for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])

        tend = my_second();
        timecomp += timediff(tstart, tend);
        sumup_large_ints(1, &npleft, &ntot);
        if(ntot > 0)
        {
            iter++;
            if(iter > 0 && ThisTask == 0) {if(iter > 10) printf("DMB: ngb iteration %d: need to repeat for %d%09d particles.\n", iter, (int) (ntot / 1000000000), (int) (ntot % 1000000000));}
            if(iter > MAXITER) {printf("DMB: failed to converge in neighbour iteration in dmb_calc()\n"); fflush(stdout); endrun(1155);}
        }
    }
    while(ntot > 0);
    
    /* iteration is done - de-malloc everything now */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    myfree(Right); myfree(Left);

    /* mark as active again */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(P[i].TimeBin < 0) {P[i].TimeBin = -P[i].TimeBin - 1;}
    }
    
    /* now that we are DONE iterating to find hsml, we can do the REAL final operations on the results */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if(dmb_isactive(i))
        {
            if (P[i].DMB_NumNgb > 0) {
                int k; for (k = 0; k < 3; k++) { P[i].DMB_V[k] /= P[i].DMB_NumNgb; }
                P[i].DMB_Density /= P[i].DMB_NumNgb;
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
    
    /* collect some timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp; CPU_Step[CPU_AGSDENSWAIT] += timewait;
    CPU_Step[CPU_AGSDENSCOMM] += timecomm; CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */

#endif
