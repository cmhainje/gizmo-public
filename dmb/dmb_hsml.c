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
    } else {
        for (k = 0; k < 3; ++k) {
            P[i].DMB_kick[k] += out->kick[k];
        }

        P[i].DMB_probtotal += out->prob_total;
        if (P[i].DMB_probtotal > 0.1) {
            printf("warning: DM particle %d had total scattering probability %e\n", P[i].ID, P[i].DMB_probtotal);
        }
    }

    MIN_ADD(P[i].DMB_dtime, out->dtime_dmb, mode);
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
    int startnode, numngb_inbox, listindex = 0, j, k, n;
    struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    if(mode == 0) {INPUTFUNCTION_NAME(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    if(local.Mass <= 0 || local.AGS_Hsml <= 0) return 0;
    int bitflag = dmb_BITFLAG(local.Type); // determine allowed particle types for search for adaptive gravitational softening terms

    /* now set particle-i centric quantities so we don't do it inside the loop */
    out.dtime_dmb = local.dtime_dmb;
    double hsml_i = local.AGS_Hsml;
    double hinv_i = 1.0 / hsml_i;
    double hinv3_i = CUBE(hinv_i);
    double hsml_j, hinv_j, hinv3_j;

    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            // search out to double the AGS_Hsml to ensure we see particles with the same AGS_Hsml as us
            double search_len = 2.001 * local.AGS_Hsml;

            numngb_inbox = ngb_treefind_pairs_threads_targeted(local.Pos, search_len, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, bitflag);
            if(numngb_inbox < 0) {return -2;} /* no neighbors! */
            for(n = 0; n < numngb_inbox; n++) /* neighbor loop */
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                hsml_j = get_hsml(j);

                if((P[j].Mass <= 0) || (hsml_j <= 0)) continue; /* make sure neighbor is valid */

                /**
                 * for all valid overlapping pairs (i, j), we are guaranteed to see it when hsml_i >= hsml_j, but we may miss it when hsml_i < hsml_j.
                 * to avoid double counting and ensure none are missed, process only the pairs where hsml_i >= hsml_j.
                 * when hsml_i == hsml_j, avoid double counting based on IDs by only processing when Type_i > Type_j.
                 * note: these considerations about double counting don't apply if j is asleep.
                 */
                if (TimeBinActive[P[j].TimeBin]) {
                    if (hsml_i < hsml_j) continue;
                    else if ((hsml_i == hsml_j) && (local.Type <= P[j].Type)) continue;
                }

                /* calculate position relative to target */
                double dx[3]; for (k = 0; k < 3; ++k) { dx[k] = local.Pos[k] - P[j].Pos[k]; }
                NEAREST_XYZ(dx[0], dx[1], dx[2], 1); // handle periodic box
                double r = sqrt(SQUARE(dx[0]) + SQUARE(dx[1]) + SQUARE(dx[2]));
                if (r > hsml_i + hsml_j) {
                    continue;
                }
                hinv_j = 1.0 / hsml_j;
                hinv3_j = CUBE(hinv_j);

                double m_chi, M_chi, v_chi[3], m_B, M_B, V_B[3], T_B;
                m_chi = All.DMB_DarkMatterMass;

                // load the relevant stuff, ensure it's all in physical CGS
                double vel_cgs = UNIT_VEL_IN_CGS / All.cf_atime;
                if (local.Type == 0) {
                    m_B = local.MolecularWeight;
                    M_B = local.Mass * UNIT_MASS_IN_CGS;
                    for (k = 0; k < 3; ++k) V_B[k] = local.Vel[k] * vel_cgs;
                    T_B = local.Temperature;

                    M_chi = P[j].Mass * UNIT_MASS_IN_CGS;
                    for (k = 0; k < 3; ++k) v_chi[k] = P[j].Vel[k] * vel_cgs;
                } else {
                    M_chi = local.Mass * UNIT_MASS_IN_CGS;
                    for (k = 0; k < 3; ++k) v_chi[k] = local.Vel[k] * vel_cgs;

                    m_B = SphP[j].DMB_MolecularWeight;
                    M_B = P[j].Mass * UNIT_MASS_IN_CGS;
                    for (k = 0; k < 3; ++k) V_B[k] = SphP[j].VelPred[k] * vel_cgs;
                    T_B = SphP[j].DMB_Temperature;
                }

                // hubble flow correction to v_chi - V_B (see ags_hsml.c:907)
                if(All.ComovingIntegrationOn) {
                    double hubble_corr = All.cf_hubble_a / All.cf_a2inv * vel_cgs;
                    double hubble_sign = (local.Type == 1) ? 1.0 : -1.0;
                    for (k = 0; k < 3; ++k) { v_chi[k] += hubble_sign * hubble_corr * dx[k]; }
                }

                double g_ij = dmb_overlap_lookup(r / hsml_i, hsml_j / hsml_i) * hinv3_j * All.cf_a3inv / CUBE(UNIT_LENGTH_IN_CGS);
                if (g_ij == 0) {
                    continue;
                }

                // handle effects on gas
                double dv = DIFF_MAG(v_chi, V_B);
                double disp_B = sqrt(T_B / m_B);
                double scrA, scrB; dmb_script_AB(dv, disp_B, &scrA, &scrB);
                double accel_coeff = (
                    m_chi / (m_chi + m_B)
                    * M_chi / m_chi
                    * g_ij
                    * scrA
                    / (UNIT_VEL_IN_CGS / UNIT_TIME_IN_CGS)
                    / All.cf_a2inv
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

                double _dt = local.Type == 0 ? local.dtime : GET_PARTICLE_TIMESTEP_IN_PHYSICAL(j);
                double _u  = 1.5 * T_B / m_B;
                if (heat_rate * _dt > _u) {
                    printf(
                        "DMB warning: large heat rate. Inputs:\n"
                        "  m_chi: %e,"
                        "  m_B: %e,"
                        "  M_chi: %e,"
                        "  M_B: %e,"
                        "  v_chi: [%e, %e, %e],"
                        "  V_B: [%e, %e, %e],"
                        "  T_B: %e,"
                        "  dv: %e,"
                        "  disp_B: %e,"
                        "  scrA: %e,"
                        "  scrB: %e,"
                        "  g_ij: %e,"
                        "  accel_coeff: %e,"
                        "  heat_rate: %e"
                        "\n",
                        m_chi,
                        m_B,
                        M_chi,
                        M_B,
                        v_chi[0], v_chi[1], v_chi[2],
                        V_B[0], V_B[1], V_B[2],
                        T_B,
                        dv,
                        disp_B,
                        scrA,
                        scrB,
                        g_ij,
                        accel_coeff,
                        heat_rate
                    );
                }

                if (local.Type == 0) {
                    for (k = 0; k < 3; ++k) { out.accel[k] += accel_coeff * (v_chi[k] - V_B[k]); }
                    out.heatrate += heat_rate;
                } else if (TimeBinActive[P[j].TimeBin]) {
                    for (k = 0; k < 3; ++k) {
                        #pragma omp atomic
                        SphP[j].DMB_Accel[k] += accel_coeff * (v_chi[k] - V_B[k]);
                    }
                    #pragma omp atomic
                    SphP[j].DMB_DtInternalEnergy += heat_rate;
                }

                // handle effects on DM

                // don't scatter inactive DM particle neighbors
                // (their scatterings were considered when they were last active)
                if (local.Type == 0 && !TimeBinActive[P[j].TimeBin]) {
                    continue;
                }

                // draw a random velocity from the local baryon MB distribution
                double v_sample[3];
                for (k = 0; k < 3; ++k) { v_sample[k] = V_B[k] + gsl_ran_gaussian(random_generator, disp_B); }
                dv = DIFF_MAG(v_chi, v_sample);

                // calculate probability to scatter
                double dt = local.Type == 1 ? local.dtime : GET_PARTICLE_TIMESTEP_IN_PHYSICAL(j); // always use DM's dt
                double prob = (M_B / m_B) * g_ij * dmb_cross_section(dv) * dv * (dt * UNIT_TIME_IN_CGS);

                double mratio = m_B / (m_chi + m_B);
                double prob_limit = DMIN(0.1, 1e-4 / (mratio * mratio));

                /* timestep condition not being met as desired, warn code to lower timestep next turn */
                if (prob > prob_limit) {
                    double new_dt = dt * (prob_limit / prob);
                    if (local.Type == 1) {
                        out.dtime_dmb = DMIN(out.dtime_dmb, new_dt);
                    } else {
                        #pragma omp critical
                        { if (new_dt < P[j].DMB_dtime) { P[j].DMB_dtime = new_dt; } }
                    }
                    // double reduction = prob_limit / prob;
                    // out.dtime_dmb = DMIN(out.dtime_dmb, local.dtime * reduction);
                    // double new_dtime_j = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(j) * reduction;
                    // #pragma omp critical (dmb_dtime_min)
                    // { if (new_dtime_j < P[j].DMB_dtime) { P[j].DMB_dtime = new_dtime_j; } }
                }

                if (local.Type == 1) {
                    out.prob_total += prob;
                } else {
                    #pragma omp atomic
                    P[j].DMB_probtotal += prob;
                }


                // roll a random number and apply scatter
                if (gsl_rng_uniform(random_generator) < prob) {
                    double Pj_dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(j);
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
                            out.kick[k] += kick[k] * All.cf_atime / UNIT_VEL_IN_CGS;
                        } else {
                            #pragma omp atomic
                            P[j].DMB_kick[k] += kick[k] * All.cf_atime / UNIT_VEL_IN_CGS;
                        }
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
        int k;
        if (P[i].Type == 0) {
            for (k = 0; k < 3; ++k) {SphP[i].DMB_Accel[k] = 0.;}
            SphP[i].DMB_DtInternalEnergy = 0.;
        } else if (P[i].Type == 1) {
            P[i].DMB_dtime = 10. * GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
            for (k = 0; k < 3; ++k) { P[i].DMB_kick[k] = 0.; }
            P[i].DMB_probtotal = 0.;
        }
    }

    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        P[i].DMB_dtime = 10. * GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
        if (P[i].Type == 0) {
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
        if (P[i].Type == 1) {
            int k; for (k = 0; k < 3; ++k) { P[i].Vel[k] += P[i].DMB_kick[k]; }
        }
    }
    
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {
        if (P[i].Type == 0) {
            int k; for (k = 0; k < 3; ++k) { P[i].GravAccel[k] += SphP[i].DMB_Accel[k]; }
            SphP[i].DtInternalEnergy += SphP[i].DMB_DtInternalEnergy;
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
