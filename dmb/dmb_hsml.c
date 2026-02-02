#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! \file dmb_hsml.c
 *  \brief dark matter-baryon interaction calculations scatterings
 *
 *  This file contains a loop modeled on the gas density computation which
 *  determines softening lengths (and appropriate correction terms)
 *  for the nearest DM particles about a gas cell and vice versa. It then
 *  computes the momentum and heat transfer rates between DM and gas.
 *  Written by Connor Hainje (connor.hainje@nyu.edu) 2023-2026, based on the AGSForce block in ags_hsml.c.
 */

#ifdef DM_DMB

#define CORE_FUNCTION_NAME dmb_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define CONDITIONFUNCTION_FOR_EVALUATION if(dmb_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

struct kernel_dmb
{
    double dp[3], dv[3], r, wk_i, wk_j, dwk_i, dwk_j, h_i, hinv_i, hinv3_i, hinv4_i, h_j, hinv_j, hinv3_j, hinv4_j;
};

/* structure for variables needed in evaluation sub-routines which must be passed from particles (sent to other processors) */
struct INPUT_STRUCT_NAME
{
    double Mass;
    double AGS_Hsml;
    double Pos[3];
    double Vel[3];
    double Temperature;
    double MolecularWeight;
    int NodeList[NODELISTLENGTH];
    int Type;
    MyIDType ID;
    double dtime;
    double dtime_dmb;
}
*DATAIN_NAME, *DATAGET_NAME;

static inline double get_hsml(MyIDType i) {
    if (P[i].Type == 0) {
        return PPP[i].Hsml;
    } else {
#ifdef AGS_HSML_CALCULATION_IS_ACTIVE
        return PPP[i].AGS_Hsml;
#else
        return All.ForceSoftening[P[i].Type];
#endif // AGS_HSML_CALCULATION_IS_ACTIVE
    }
}

/* routine to pass particle information to the actual evaluation sub-routines */
static inline void INPUTFUNCTION_NAME(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    in->Mass = PPP[i].Mass;

    in->AGS_Hsml = get_hsml(i);

    if (P[i].Type == 0) {
        in->Temperature = SphP[i].DMB_Temperature;
        in->MolecularWeight = SphP[i].DMB_MolecularWeight;
    } else {
        in->Temperature = 0.;
        in->MolecularWeight = 0.;
    }

    in->Type = P[i].Type;
    in->ID = P[i].ID;
    in->dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
    in->dtime_dmb = P[i].DMB_dtime;
    int k; for(k = 0; k < 3; ++k) {
        in->Pos[k] = P[i].Pos[k];
        in->Vel[k] = P[i].Vel[k];
    }
}


/* structure for variables which must be returned -from- the evaluation sub-routines */
struct OUTPUT_STRUCT_NAME
{
    // gas outputs
    double accel[3];
    double heatrate;

    // DM outputs
    double kick[3];
    double dtime_dmb;
    double prob_total;
    
    int numngb;
    int numscat;
}
*DATARESULT_NAME, *DATAOUT_NAME;

#define ASSIGN_ADD_PRESET(x,y,mode) (mode == 0 ? (x=y) : (x+=y))
#define MINMAX_CHECK(x,xmin,xmax) ((x<xmin)?(xmin=x):((x>xmax)?(xmax=x):(1)))
#define MAX_ADD(x,y,mode) ((y > x) ? (x = y) : (1)) // simpler definition now used
#define MIN_ADD(x,y,mode) ((y < x) ? (x = y) : (1))

static inline void OUTPUTFUNCTION_NAME(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    int k; 
    if (P[i].Type == 0) {
        for (k = 0; k < 3; ++k) {SphP[i].DMB_Accel[k] += out->accel[k];}
        SphP[i].DMB_DtInternalEnergy += out->heatrate;

        for (k = 0; k < 3; ++k) {P[i].GravAccel[k] += SphP[i].DMB_Accel[k];}
        SphP[i].DtInternalEnergy += SphP[i].DMB_DtInternalEnergy;
    } else {
        for (k = 0; k < 3; ++k) {
            P[i].DMB_kick[k] += out->kick[k];
            // P[i].Vel[k] += P[i].DMB_kick[k];
        }

        // for (k = 0; k < 3; ++k) {
        //     P[i].Vel[k] += out->kick[k];
        // }

        P[i].DMB_probtotal += out->prob_total;
        if (P[i].DMB_probtotal > 0.2) {
            printf("warning: DM particle %d had total scattering probability %e\n", P[i].ID, P[i].DMB_probtotal);
        }

        P[i].DMB_NumNeighbors += out->numngb;
        P[i].DMB_NumScatters += out->numscat;
        // if (P[i].DMB_NumScatters > 0) {
        //     printf("  %d scattered with %d / %d neighbors\n", P[i].ID, P[i].DMB_NumScatters, P[i].DMB_NumNeighbors);
        // }


        // for (k = 0; k < 3; ++k) {P[i].Vel[k] += out->kick[k];}
        MIN_ADD(P[i].DMB_dtime, out->dtime_dmb, mode);
    }
}


/* routine to determine if we need to apply the additional DM-b calculation[s] */
int dmb_isactive(int i);
int dmb_isactive(int i)
{
    if (P[i].TimeBin < 0) return 0; /* check our 'marker' for particles which have finished iterating to an Hsml solution (if they have, dont do them again) */
    if (P[i].Mass < 0) return 0;
    if (P[i].Type <= 1) return 1; /* accept gas and DM */
    return 0;
}

int dmb_BITFLAG(short int particle_type_primary)
{
    if (particle_type_primary == 0) {return 2;} /* 2^1: gas sees DM */
    if (particle_type_primary == 1) {return 1;} /* 2^0: DM sees gas */
    return 0;
}

/* generates a random unit-length 3vector */
void random_unit_vector(double out[3]) {
    double cos_theta = 2.0 * gsl_rng_uniform(random_generator) - 1.0;
    double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
    double phi = gsl_rng_uniform(random_generator) * 2.0 * M_PI;
    out[0] = sin_theta * cos(phi);
    out[1] = sin_theta * sin(phi);
    out[2] = cos_theta;
}


#define SQUARE(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define DIFF_MAG(v1, v2) (sqrt(SQUARE(v1[0]-v2[0]) + SQUARE(v1[1]-v2[1]) + SQUARE(v1[2]-v2[2])))

/*!   -- this subroutine writes to shared memory [updating the neighbor values]: need to protect these writes for openmp below. none of the modified values are read, so only the write block is protected. note the writes can occur in the called code-blocks, so need to make sure they are followed so everything can be carefully constructed */
int dmb_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    /* zero memory and import data for local target */
    int startnode, numngb_inbox, listindex = 0, j, k, n; double r2, u_i, u_j;
    struct kernel_dmb kernel; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME)); memset(&kernel, 0, sizeof(struct kernel_dmb));
    if(mode == 0) {INPUTFUNCTION_NAME(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    if(local.Mass <= 0 || local.AGS_Hsml <= 0) return 0;
    /* now set particle-i centric quantities so we don't do it inside the loop */
    kernel.h_i = local.AGS_Hsml; kernel_hinv(kernel.h_i, &kernel.hinv_i, &kernel.hinv3_i, &kernel.hinv4_i);
    int bitflag = dmb_BITFLAG(local.Type); // determine allowed particle types for search for adaptive gravitational softening terms
    out.dtime_dmb = local.dtime_dmb;

    // printf("type: %d, hsml: %e\n", local.Type, local.AGS_Hsml);


    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            // search out to at least double the AGS_Hsml to ensure we see particles with the same AGS_Hsml as us
            double search_len = 2.5 * local.AGS_Hsml;

            numngb_inbox = ngb_treefind_pairs_threads_targeted(local.Pos, search_len, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, bitflag);
            if(numngb_inbox < 0) {return -2;} /* no neighbors! */
            for(n = 0; n < numngb_inbox; n++) /* neighbor loop */
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                double Pj_Hsml = get_hsml(j);

                if((P[j].Mass <= 0) || (Pj_Hsml <= 0)) continue; /* make sure neighbor is valid */

                /**
                 * for all valid overlapping pairs (i, j), we are guaranteed to see it when hsml_i >= hsml_j, but we may miss it when hsml_i < hsml_j.
                 * to avoid double counting and ensure none are missed, process only the pairs where hsml_i >= hsml_j.
                 * when hsml_i == hsml_j, avoid double counting based on IDs by only processing when Type_i > Type_j.
                 * note: these considerations about double counting don't apply if j is asleep.
                 */

                if (TimeBinActive[P[j].TimeBin]) {
                    if (local.AGS_Hsml < Pj_Hsml) continue;
                    else if ((local.AGS_Hsml == Pj_Hsml) && (local.Type <= P[j].Type)) continue;
                }

                /* calculate position relative to target */
                for (k = 0; k < 3; ++k) { kernel.dp[k] = local.Pos[k] - P[j].Pos[k]; }
                NEAREST_XYZ(kernel.dp[0], kernel.dp[1], kernel.dp[2], 1); /*  now find the closest image in the given box size  */
                r2 = kernel.dp[0]*kernel.dp[0] + kernel.dp[1]*kernel.dp[1] + kernel.dp[2]*kernel.dp[2];
                if (r2 <= 0) continue;
                kernel.r = sqrt(r2);
                kernel.h_j = Pj_Hsml;
                if (kernel.r > kernel.h_i + kernel.h_j) continue;
                /* calculate kernel quantities needed below */
                kernel_hinv(kernel.h_j, &kernel.hinv_j, &kernel.hinv3_j, &kernel.hinv4_j);
                u_i = kernel.r * kernel.hinv_i; u_j = kernel.r * kernel.hinv_j;
                if(u_i < 1) {kernel_main(u_i, kernel.hinv3_i, kernel.hinv4_i, &kernel.wk_i, &kernel.dwk_i, 0);} else {kernel.wk_i=kernel.dwk_i=0;}
                if(u_j < 1) {kernel_main(u_j, kernel.hinv3_j, kernel.hinv4_j, &kernel.wk_j, &kernel.dwk_j, 0);} else {kernel.wk_j=kernel.dwk_j=0;}
                for(k = 0; k < 3; ++k) {
                    kernel.dv[k] = local.Vel[k] - P[j].Vel[k];
                    if (All.ComovingIntegrationOn) { kernel.dv[k] += All.cf_hubble_a * kernel.dp[k] / All.cf_a2inv; }
                }

                double m_chi, M_chi, v_chi[3], m_B, M_B, V_B[3], T_B;
                m_chi = All.DMB_DarkMatterMass;

                // load the relevant stuff, ensure it's all in CGS
                if (local.Type == 0) {
                    m_B = local.MolecularWeight;
                    M_B = local.Mass * UNIT_MASS_IN_CGS;
                    for (k = 0; k < 3; ++k) V_B[k] = local.Vel[k] * UNIT_VEL_IN_CGS;
                    T_B = local.Temperature;

                    M_chi = P[j].Mass * UNIT_MASS_IN_CGS;
                    for (k = 0; k < 3; ++k) v_chi[k] = P[j].Vel[k] * UNIT_VEL_IN_CGS;
                } else {
                    M_chi = local.Mass * UNIT_MASS_IN_CGS;
                    for (k = 0; k < 3; ++k) v_chi[k] = local.Vel[k] * UNIT_VEL_IN_CGS;

                    m_B = SphP[j].DMB_MolecularWeight;
                    M_B = P[j].Mass * UNIT_MASS_IN_CGS;
                    for (k = 0; k < 3; ++k) V_B[k] = P[j].Vel[k] * UNIT_VEL_IN_CGS;
                    T_B = SphP[j].DMB_Temperature;
                }

                double g_ij = dmb_overlap_lookup(
                    kernel.r / kernel.h_i,
                    kernel.h_j / kernel.h_i
                ) * kernel.hinv3_j / CUBE(UNIT_LENGTH_IN_CGS);

                if (g_ij == 0) {
                    continue;
                }

                // handle effects on gas
                double dv = DIFF_MAG(v_chi, V_B);
                double scrA = dmb_script_A(dv, T_B / m_B);
                double scrB = dmb_script_B(dv, T_B / m_B);
                double accel_coeff = (
                    m_chi / (m_chi + m_B)
                    * M_chi / m_chi
                    * g_ij
                    * scrA
                    / (UNIT_VEL_IN_CGS / UNIT_TIME_IN_CGS)
                );
                double heat_rate = (
                    m_chi / (m_chi + m_B)
                    * M_chi / m_chi
                    * g_ij
                    * (
                        - m_B / (m_chi + m_B) * scrB
                        + dv * dv * scrA
                    )
                    / (UNIT_SPECEGY_IN_CGS / UNIT_TIME_IN_CGS)
                );

                // if (P[j].ID == 32768) {
                //     __builtin_trap();
                // }

                if (local.Type == 0) {
                    for (k = 0; k < 3; ++k) { 
                        out.accel[k] += accel_coeff * (v_chi[k] - V_B[k]);
                    }
                    out.heatrate += heat_rate;
                } else {
                    for (k = 0; k < 3; ++k) { 
                        #pragma omp atomic
                        SphP[j].DMB_Accel[k] += accel_coeff * (v_chi[k] - V_B[k]);
                    }
                    #pragma omp atomic
                    SphP[j].DMB_DtInternalEnergy += heat_rate;
                }


                // handle effects on DM

                // draw a random velocity from the local baryon MB distribution
                double v_sample[3];
                for (k = 0; k < 3; ++k) { v_sample[k] = V_B[k] + gsl_ran_gaussian(random_generator, sqrt(T_B / m_B)); }
                dv = DIFF_MAG(v_chi, v_sample);

                // calculate probability to scatter
                // double dt = ((local.Type == 1) ? local.dtime : GET_PARTICLE_TIMESTEP_IN_PHYSICAL(j)) * UNIT_TIME_IN_CGS;
                double prob = (M_B / m_B) * g_ij * dmb_cross_section(dv) * dv * (local.dtime * UNIT_TIME_IN_CGS);
                // double prob = (M_B / m_B) * g_ij * scrA * (local.dtime * UNIT_TIME_IN_CGS);
                if (prob > 0.2) {
                    printf("warning: large probability to scatter (IDs: %d, %d, prob: %e)\n", local.ID, P[j].ID, prob);
                    out.dtime_dmb = DMIN(out.dtime_dmb , local.dtime * (0.2 / prob)); /* timestep condition not being met as desired, warn code to lower timestep next turn */
                }

                if (local.Type == 1) {
                    out.prob_total += prob;
                    out.numngb += 1;
                } else {
                    #pragma omp atomic
                    P[j].DMB_probtotal += prob;

                    #pragma omp atomic
                    P[j].DMB_NumNeighbors += 1;
                }


                // roll a random number and apply scatter
                double Pj_dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(j);
                if (gsl_rng_uniform(random_generator) < prob) {
#ifdef WAKEUP
                    if (!(TimeBinActive[P[j].TimeBin])) {
                        if (WAKEUP * local.dtime < Pj_dtime) {
                            #pragma omp atomic write
                            PPPZ[j].wakeup = 1;
                            #pragma omp atomic write
                            NeedToWakeupParticles_local = 1;
                        }
                    }
#endif
                    double ehat[3]; random_unit_vector(ehat);

                    double kick[3];
                    for (k = 0; k < 3; ++k) {
                        kick[k] = m_B / (m_chi + m_B) * (
                            v_sample[k] - v_chi[k]
                            + dv * ehat[k]
                        );

                        if (local.Type == 1) {
                            out.kick[k] += kick[k] / UNIT_VEL_IN_CGS;
                        } else {
                            #pragma omp atomic
                            P[j].DMB_kick[k] += kick[k] / UNIT_VEL_IN_CGS;
                            // P[j].Vel[k] += kick[k] / UNIT_VEL_IN_CGS;
                        }
                    }

                    if (local.Type == 1) {
                        out.numscat += 1;
                    } else {
                        #pragma omp atomic
                        P[j].DMB_NumScatters += 1;
                    }
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
    PRINT_STATUS(" ..entering DM-b calculation");
    /* before doing any operations, need to zero the appropriate memory so we can correctly do pair-wise operations */
    int i;
    for(i = 0; i < NumPart; ++i) {
    // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        int k;
        if (P[i].Type == 1) {
            P[i].DMB_dtime = 10. * GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
            for (k = 0; k < 3; ++k) { P[i].DMB_kick[k] = 0.; }
            P[i].DMB_probtotal = 0.;
            P[i].DMB_NumScatters = 0;
            P[i].DMB_NumNeighbors = 0;
        }
    }

    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        int k;
        P[i].DMB_dtime = 10. * GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
        for (k = 0; k < 3; ++k) { P[i].DMB_kick[k] = 0.; }
        P[i].DMB_probtotal = 0.;
        P[i].DMB_NumScatters = 0;
        P[i].DMB_NumNeighbors = 0;
        
        if (P[i].Type == 0) {
            for (k = 0; k < 3; ++k) {SphP[i].DMB_Accel[k] = 0.;}
            SphP[i].DMB_DtInternalEnergy = 0.;

            double u = SphP[i].InternalEnergyPred;
            double mu, T;
            #ifdef COOLING
            double rho = SphP[i].Density * All.cf_a3inv;
            double ne=1, nh0=0, nHe0, nHepp, nhp, nHeII;
            T = ThermalProperties(u, rho, i, &mu, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp);
            #else
            mu = 1.0;
            T = u * ((2./3.) * mu * U_TO_TEMP_UNITS);
            #endif // COOLING

            SphP[i].DMB_Temperature = T * BOLTZMANN_CGS;
            SphP[i].DMB_MolecularWeight = mu * PROTONMASS_CGS;
        }
    }
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    /* do final operations on results: these are operations that can be done after the complete set of iterations */
    for(i = 0; i < NumPart; ++i) {
    // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        int k;
        if (P[i].Type == 1) {
            for (k = 0; k < 3; ++k) {
                P[i].Vel[k] += P[i].DMB_kick[k];
            }
        }
    }

    /* collect timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp; CPU_Step[CPU_AGSDENSWAIT] += timewait;
    CPU_Step[CPU_AGSDENSCOMM] += timecomm; CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);
    // __builtin_trap();
    PRINT_STATUS(" ..DM-b finished (%f sec)", timeall);
    // endrun(1234);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */

#endif // DM_DMB
