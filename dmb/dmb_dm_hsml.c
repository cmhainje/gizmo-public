#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! \file dmb_dm_hsml
 *  \brief smoothing length and velocity dispersion calculation for dark matter particles around gas and dark matter particles
 *
 *  This file contains a loop modeled on the gas density computation which
 *    determines kernel lengths for dark matter particles around a given set of gas and dark matter particles;
 *    this is used by the flag DM_DMB to estimate the local dark matter velocity dispersion around a given gas 
 *    or dark matter particle, which (in turn) is used to compute the momentum and heat exchange rates between
 *    dark matter and baryons due to scattering. The loop here needs to be called for these models (note this
 *    in general will require a different smoothing length from, say, the dm-dm force softening, or the
 *    gas kernel length for hydro, hence it requires a whole additional loop, even though the loop is
 *    functionally identical - modulo which particles are used - to the loop for adaptive gravitational softening
 *
 * This file was originally written by Qirong Zhu, for GIZMO, based on Phil Hopkins's adaptive gravitational softening
 *    routine. It has been here copied and modified by Connor Hainje.
 *
 */

#ifdef DM_DMB

#define CORE_FUNCTION_NAME disp_density_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME disp_particle2in_density    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME disp_out2particle_density  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(disp_density_isactive(i)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

/* this structure defines the variables that need to be sent -from- the 'searching' element */
static struct INPUT_STRUCT_NAME
{
    MyDouble Pos[3];
    MyFloat Hsml;
    int NodeList[NODELISTLENGTH];
}
*DATAIN_NAME, *DATAGET_NAME;

/* this subroutine assigns the values to the variables that need to be sent -from- the 'searching' element */
void disp_particle2in_density(struct INPUT_STRUCT_NAME *in, int i, int loop_iteration)
{
    int k; for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k];}
    in->Hsml = P[i].DMB_Hsml;
}


/* this structure defines the variables that need to be sent -back to- the 'searching' element */
static struct OUTPUT_STRUCT_NAME
{
    MyLongDouble DM_Ngb;
    MyLongDouble DM_Vx;
    MyLongDouble DM_Vy;
    MyLongDouble DM_Vz;
    MyLongDouble DM_VelDisp;
    MyLongDouble DM_Density;

    MyLongDouble Gas_Ngb;
    MyLongDouble Gas_Vx;
    MyLongDouble Gas_Vy;
    MyLongDouble Gas_Vz;
    MyLongDouble Gas_Density;
    MyLongDouble Gas_Temperature;
    MyLongDouble Gas_MicroMass;
}
*DATARESULT_NAME, *DATAOUT_NAME;

/* this subroutine assigns the values to the variables that need to be sent -back to- the 'searching' element */
void disp_out2particle_density(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    ASSIGN_ADD(P[i].DMB_NumNgb, out->DM_Ngb + out->Gas_Ngb, mode);

    ASSIGN_ADD(P[i].DMB_NumNgbDM, out->DM_Ngb, mode);
    ASSIGN_ADD(P[i].DMB_VDM[0], out->DM_Vx, mode);
    ASSIGN_ADD(P[i].DMB_VDM[1], out->DM_Vy, mode);
    ASSIGN_ADD(P[i].DMB_VDM[2], out->DM_Vz, mode);
    ASSIGN_ADD(P[i].DMB_VelDispDM, out->DM_VelDisp, mode);
    ASSIGN_ADD(P[i].DMB_DensityDM, out->DM_Density, mode);

    ASSIGN_ADD(P[i].DMB_NumNgbGas, out->Gas_Ngb, mode);
    ASSIGN_ADD(P[i].DMB_VGas[0], out->Gas_Vx, mode);
    ASSIGN_ADD(P[i].DMB_VGas[1], out->Gas_Vy, mode);
    ASSIGN_ADD(P[i].DMB_VGas[2], out->Gas_Vz, mode);
    ASSIGN_ADD(P[i].DMB_DensityGas, out->Gas_Density, mode);
    ASSIGN_ADD(P[i].DMB_TemperatureGas, out->Gas_Temperature, mode);
    ASSIGN_ADD(P[i].DMB_MicroparticleMassGas, out->Gas_MicroMass, mode);
}


/* routine to determine if we need to use disp_density to calculate Hsml */
int disp_density_isactive(int i);
int disp_density_isactive(int i)
{
    if(P[i].TimeBin < 0) return 0;
    if(P[i].Type > 1) return 0; // only gas and DM particles //
    if(P[i].Mass <= 0) return 0;
    return 1;
}


/*! This function represents the core of the density computation. The target particle may either be local, or reside in the communication buffer. */
/*!   -- this subroutine contains no writes to shared memory -- */
int disp_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, n; struct INPUT_STRUCT_NAME local; struct OUTPUT_STRUCT_NAME out; memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME)); /* define variables and zero memory and import data for local target*/
    if(mode == 0) {INPUTFUNCTION_NAME(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];} /* imports the data to the correct place and names */
    /* Now start the actual neighbor computation for this particle */
    if(mode == 0) {startnode = All.MaxPart; /* root node */} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;    /* open it */}
    while(startnode >= 0) {
        while(startnode >= 0) {
            numngb_inbox = ngb_treefind_variable_threads_targeted(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist, 3); // search for high-res DM and gas particles: 2^0 + 2^1 = 3
            if(numngb_inbox < 0) {return -2;}
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if(P[j].Mass <= 0) continue;

                // Gas
                if (P[j].Type == 0) {
                    out.Gas_Vx += P[j].Vel[0];
                    out.Gas_Vy += P[j].Vel[1];
                    out.Gas_Vz += P[j].Vel[2];
                    out.Gas_Density += P[j].Mass;

                    double rho_B = SphP[j].Density * All.cf_a3inv;
                    double u_B = SphP[j].InternalEnergyPred;
                    double mu=1, ne=1, nh0=0, nHe0, nHepp, nhp, nHeII; // pull various known thermal properties, prepare to extract others //
                    double T_B = ThermalProperties(u_B, rho_B, j, &mu, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties
                    out.Gas_Temperature += T_B;
                    out.Gas_MicroMass += mu;

                    out.Gas_Ngb++;
                }

                // Dark matter
                else if (P[j].Type == 1) {
                    out.DM_Vx += P[j].Vel[0];
                    out.DM_Vy += P[j].Vel[1];
                    out.DM_Vz += P[j].Vel[2];
                    out.DM_VelDisp += (P[j].Vel[0] * P[j].Vel[0] + P[j].Vel[1] * P[j].Vel[1] + P[j].Vel[2] * P[j].Vel[2]);
                    out.DM_Density += P[j].Mass;
                    out.DM_Ngb++;
                }
            } // numngb_inbox loop
        } // while(startnode)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode; /* open it */}}} /* continue to open leaves if needed */
    }
    if(mode == 0) {OUTPUTFUNCTION_NAME(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;} /* collects the result at the right place */
    return 0;
}


void disp_density(void)
{
    /* initialize variables used below, in particlar the structures we need to call throughout the iteration */
    CPU_Step[CPU_MISC] += measure_time(); double t00_truestart = my_second();
    MyFloat *Left, *Right; double desnumngb=64, desnumngbdev=48; long long ntot; 
    int i, npleft, iter=0, redo_particle;
    Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
    Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));
    /* initialize anything we need to about the active particles before their loop */
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) {if(disp_density_isactive(i)) {P[i].DMB_NumNgb = 0; Left[i] = Right[i] = 0;}}
    
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
            if(disp_density_isactive(i))
            {
                redo_particle = 0; /* now check whether we have enough neighbours, and are below the maximum search radius */
                double maxsoft = DMIN(All.MaxHsml, 10.0*PPP[i].Hsml);
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
                        PRINT_WARNING("DM disp: i=%d task=%d ID=%llu Type=%d Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
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
                            char buf[1000]; sprintf(buf, "DM disp: Right[i] == 0 && Left[i] == 0 && P[i].DMB_Hsml=%g\n", P[i].DMB_Hsml); terminate(buf);
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
            } //  if(disp_density_isactive(i))
        } // for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])

        tend = my_second();
        timecomp += timediff(tstart, tend);
        sumup_large_ints(1, &npleft, &ntot);
        if(ntot > 0)
        {
            iter++;
            if(iter > 0 && ThisTask == 0) {if(iter > 10) printf("DM disp: ngb iteration %d: need to repeat for %d%09d particles.\n", iter, (int) (ntot / 1000000000), (int) (ntot % 1000000000));}
            if(iter > MAXITER) {printf("DM disp: failed to converge in neighbour iteration in disp_density()\n"); fflush(stdout); endrun(1155);}
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
        int k;
        if(disp_density_isactive(i))
        {
            // Gas
            if (P[i].DMB_NumNgbGas > 0) {
                for (k = 0; k < 3; k++) { P[i].DMB_VGas[k] /= P[i].DMB_NumNgbGas; }
                P[i].DMB_TemperatureGas /= P[i].DMB_NumNgbGas;
                P[i].DMB_MicroparticleMassGas /= P[i].DMB_NumNgbGas;
                P[i].DMB_DensityGas /= P[i].DMB_Hsml * P[i].DMB_Hsml * P[i].DMB_Hsml;
            } else {
                if (P[i].DMB_TemperatureGas < 0 || isnan(P[i].DMB_TemperatureGas)) {
                    P[i].DMB_TemperatureGas = 0;
                }
                if (P[i].DMB_MicroparticleMassGas < 0 || isnan(P[i].DMB_MicroparticleMassGas)) {
                    P[i].DMB_MicroparticleMassGas = 0;
                }
                if (P[i].DMB_DensityGas < 0 || isnan(P[i].DMB_DensityDM)) {
                    P[i].DMB_DensityGas = 0;
                }
            }

            // Dark matter
            if (P[i].DMB_NumNgbDM > 0) {
                for (k = 0; k < 3; k++) { P[i].DMB_VDM[k] /= P[i].DMB_NumNgbDM; }
                P[i].DMB_VelDispDM /= P[i].DMB_NumNgbDM;
                P[i].DMB_VelDispDM = (1./All.cf_atime) * sqrt(P[i].DMB_VelDispDM - P[i].DMB_VDM[0] * P[i].DMB_VDM[0] - P[i].DMB_VDM[1] * P[i].DMB_VDM[1] - P[i].DMB_VDM[2] * P[i].DMB_VDM[2]) / 1.732;  // 1d velocity dispersion
                P[i].DMB_DensityDM /= P[i].DMB_Hsml * P[i].DMB_Hsml * P[i].DMB_Hsml;
            } else {
                if (P[i].DMB_VelDispDM < 0 || isnan(P[i].DMB_VelDispDM)) {
                    P[i].DMB_VelDispDM = sqrt(P[i].Vel[0]*P[i].Vel[0]+P[i].Vel[1]*P[i].Vel[1]+P[i].Vel[2]*P[i].Vel[2])/All.cf_atime;
                }
                if (P[i].DMB_DensityDM < 0 || isnan(P[i].DMB_DensityDM)) {
                    P[i].DMB_DensityDM = 0;
                }
            }

            compute_exch_rates(i, P[i].DMB_MomExch, &P[i].DMB_HeatExch);

            double dtime = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i); /*  the actual time-step */

            // current particle is gas
            if (P[i].Type == 0) {
                // momentum kick
                double kick_coeff = dtime / P[i].DMB_DensityGas * All.cf_atime;
                for (k = 0; k < 3; k++) {
                    SphP[i].VelPred[k] += P[i].DMB_MomExch[k] * kick_coeff;
                    P[i].Vel[k] += P[i].DMB_MomExch[k] * kick_coeff;
                }

                // heat kick
                double volume = P[i].DMB_Hsml * P[i].DMB_Hsml * P[i].DMB_Hsml; /* is there a better estimate? */
                SphP[i].InternalEnergy += P[i].DMB_HeatExch * dtime * volume;
                SphP[i].InternalEnergyPred += P[i].DMB_HeatExch * dtime * volume;
            }

            // current particle is DM
            else if (P[i].Type == 1) {
                // momentum kick only
                double kick_coeff = dtime / P[i].DMB_DensityDM * All.cf_atime;
                for (k = 0; k < 3; k++) { P[i].Vel[k] += P[i].DMB_MomExch[k] * kick_coeff; }
            }
        }
    }
    
    /* collect some timing information */
    double t1; t1 = WallclockTime = my_second(); timeall = timediff(t00_truestart, t1);
    CPU_Step[CPU_AGSDENSCOMPUTE] += timecomp; CPU_Step[CPU_AGSDENSWAIT] += timewait;
    CPU_Step[CPU_AGSDENSCOMM] += timecomm; CPU_Step[CPU_AGSDENSMISC] += timeall - (timecomp + timewait + timecomm);
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */


#endif




