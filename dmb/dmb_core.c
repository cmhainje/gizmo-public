#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_randist.h>

#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! \file dmb_core.c
 *  \brief Fuctions needed for the calculations of DM-baryon interactions
 *
 *  This file contains the functions and routines necesary for the calculation of
 *  the effect of interactions between dark matter and baryons.
 *  Written by Connor Hainje (connor.hainje@nyu.edu) 2023-2026.
 */

#ifdef DM_DMB

#define GSLWORKSIZE 100000

double HYPERG_ASYMP_FACTOR = 100; /*!< use asymptotic limit when inputs exceed this factor */

/**
 * Computes the value of the cross section at a given relative velocity `v`.
 * Units:
 *   v:      physical, cgs [cm/s]
 *   return: physical, cgs [cm^2]
*/
double dmb_cross_section(double v) {
  return All.DMB_InteractionCrossSection * pow(v / C_LIGHT_CGS, All.DMB_InteractionPowerScale);
}

/**
 * Implements the function \mathcal{A}.
 * Units:
 *   v:      physical, cgs [cm/s]
 *   disp:   physical, cgs [cm^2/s^2]
 *   return: physical, cgs [cm^3/s]
 */
double dmb_script_A(double w, double disp) {
  if (w * w > HYPERG_ASYMP_FACTOR * disp) { return dmb_cross_section(w) * w; }

  int n = All.DMB_InteractionPowerScale;
  double sqrt_disp = sqrt(disp);
  double out = (
    sqrt(pow(2.0, 5.0 + n) / M_PI) / 3.0
    * gsl_sf_gamma(3.0 + 0.5 * n)
    * dmb_cross_section(sqrt_disp) * sqrt_disp
    * gsl_sf_hyperg_1F1(-0.5 * (n + 1), 2.5, -0.5 * w * w / disp)
  );
  if (isnan(out)) {
    printf("ERROR: script_A returning NaN; inputs were w=%f, disp=%f\n", w, disp);
    endrun(9999);
  }
  return out;
}

/**
 * Implements the function \mathcal{B}.
 * Units:
 *   v:      physical, cgs [cm/s]
 *   disp:   physical, cgs [cm^2/s^2]
 *   return: physical, cgs [cm^3/s]
 */
double dmb_script_B(double w, double disp) {
  if (w * w > HYPERG_ASYMP_FACTOR * disp) { return dmb_cross_section(w) * w * w * w; }

  int n = All.DMB_InteractionPowerScale;
  double sqrt_disp = sqrt(disp);
  double out = (
    sqrt(pow(2.0, 5.0 + n) / M_PI)
    * gsl_sf_gamma(3.0 + 0.5 * n)
    * dmb_cross_section(sqrt_disp) * sqrt_disp * sqrt_disp * sqrt_disp
    * gsl_sf_hyperg_1F1(-0.5 * (n + 3), 1.5, -0.5 * w * w / disp)
  );
  if (isnan(out)) {
    printf("ERROR: script_B returning NaN; inputs were w=%f, disp=%f\n", w, disp);
    endrun(9999);
  }
  return out;
}

void dmb_script_AB(double w, double disp, double *scrA, double *scrB) {
  int n = All.DMB_InteractionPowerScale;
  double sigma = dmb_cross_section(disp);

  if (w * w > HYPERG_ASYMP_FACTOR * disp * disp) {
    *scrA = sigma * w;
    *scrB = sigma * w * w * w;
    return;
  }

  *scrA = (
    All.DMB_ScriptCoeff / 3.0 * sigma * disp
    * gsl_sf_hyperg_1F1(-0.5 * (n + 1), 2.5, -0.5 * w * w / (disp * disp))
  );
  if (isnan(*scrA)) {
    printf("ERROR: script_A is NaN; inputs were w=%f, disp=%f\n", w, disp);
    endrun(9999);
  }

  *scrB = (
    All.DMB_ScriptCoeff * sigma * disp * disp * disp
    * gsl_sf_hyperg_1F1(-0.5 * (n + 3), 1.5, -0.5 * w * w / (disp * disp))
  );
  if (isnan(*scrB)) {
    printf("ERROR: script_B is NaN; inputs were w=%f, disp=%f\n", w, disp);
    endrun(9999);
  }
}

/*! initialize some variables */
void dmb_init() {
  int i, k;
  for (i = 0; i < NumPart; ++i) {
    if (P[i].Type == 0) {
      for (k = 0; k < 3; ++k) { SphP[i].DMB_Accel[k] = 0.; }
      SphP[i].DMB_DtInternalEnergy = 0.;
      SphP[i].DMB_MolecularWeight = 0.;
      SphP[i].DMB_Temperature = 0.;
    }
    P[i].DMB_dtime = 0;
  }

  All.DMB_ScriptCoeff = sqrt(pow(2.0, 5.0 + All.DMB_InteractionPowerScale) / M_PI) * gsl_sf_gamma(3.0 + 0.5 * All.DMB_InteractionPowerScale);
}



/* *** HANDLING THE OVERLAP INTEGRAL *** */

/**
 * calculate the index of the given h_ratio
 * (h_ratios are logarithmically spaced from 1 to 0.01)
 */
int h_ratio_index(double h_ratio) {
  double log_hr = log10(h_ratio > 1 ? 1 / h_ratio : h_ratio);
  double raw_index = -0.5 * log_hr * DMB_OVERLAP_NUM_H_RATIO;
  return DMIN(floor(raw_index), DMB_OVERLAP_NUM_H_RATIO - 2);
}

/**
 * calculate the index of the given delta
 * (deltas are linearly spaced from 0 to 2)
 */
int delta_index(double delta) {
  double raw_index = 0.5 * delta * DMB_OVERLAP_NUM_DELTA;
  return DMIN(floor(raw_index), DMB_OVERLAP_NUM_DELTA - 2);
}


#define CLAMP(x, lo, hi) DMAX(DMIN(x, hi), lo)

/**
 * linearly interpolate to coordinate t given values a at t=0 and b at t=1
 */
double lerp(double t, double a, double b) {
  return a + CLAMP(t, 0., 1.) * (b - a);
}

/**
 * bilinearly interpolate to coordinate (x, y) given values at the four
 * quadrants of the unit square (where `Q{x}{y}` is at (x, y))
 */
double blerp(double x, double y, double Q00, double Q10, double Q01, double Q11) {
  return lerp(y, lerp(x, Q00, Q10), lerp(x, Q01, Q11));
}

/**
 * look up the value of the overlap integral in the precomputed table
 */
double dmb_overlap_lookup(double delta, double h_ratio) {
  if (h_ratio <= 0) {
    printf("ERROR: overlap integral called with h_ratio <= 0.\n");
    endrun(9999);
  }
  if (delta < 0) {
    printf("ERROR: overlap integral called with delta < 0.\n");
    endrun(9999);
  }
  if (delta >= 1 + h_ratio) {
    return 0.0;
  }
  if (h_ratio > 1) {
    return dmb_overlap_lookup(delta / h_ratio, 1 / h_ratio) * h_ratio * h_ratio * h_ratio;
  }

  double hr_index = -0.5 * log10(h_ratio) * (DMB_OVERLAP_NUM_H_RATIO - 1);
  double dl_index = 0.5 * delta * (DMB_OVERLAP_NUM_DELTA - 1);

  int i = DMIN(floor(h_ratio_index(h_ratio)), DMB_OVERLAP_NUM_H_RATIO - 2);
  int j = DMIN(floor(delta_index(delta)), DMB_OVERLAP_NUM_DELTA - 2);

  return blerp(
    hr_index - i,
    dl_index - j,
    DMB_OverlapTable[ i * DMB_OVERLAP_NUM_DELTA + j ],
    DMB_OverlapTable[ (i+1) * DMB_OVERLAP_NUM_DELTA + j ],
    DMB_OverlapTable[ i * DMB_OVERLAP_NUM_DELTA + (j+1) ],
    DMB_OverlapTable[ (i+1) * DMB_OVERLAP_NUM_DELTA + (j+1) ]
  );
}

double overlap_c_integrand(double c, void* params) {
  // unpack params
  double h_ratio, delta, s;
  h_ratio = *(double *) params;
  delta = *(double *) (params + sizeof(double));
  s = *(double *) (params + 2 * sizeof(double));

  // argument to kernel
  double f = sqrt(s * s + delta * delta - 2 * s * delta * c) / h_ratio;

  // evaluate kernel
  double wk; kernel_main(f, 1, 1, &wk, &wk, -1);

  // return the integrand
  return wk;
}

double overlap_s_integrand(double s, void* params) {
  // unpack params
  double result, abserr, h_ratio, delta, newparams[3];
  h_ratio = *(double *) params;
  delta = *(double *) (params + sizeof(double));

  // set up integration workspace
  gsl_function F;
  gsl_integration_workspace *workspace;
  workspace = gsl_integration_workspace_alloc(GSLWORKSIZE);

  // register integrand and parameters
  newparams[0] = h_ratio;
  newparams[1] = delta;
  newparams[2] = s;
  F.function = &overlap_c_integrand;
  F.params = newparams;

  // perform the integral
  gsl_integration_qag(&F, -1.0, +1.0, 0, 1e-8, GSLWORKSIZE, GSL_INTEG_GAUSS41, workspace, &result, &abserr);
  gsl_integration_workspace_free(workspace);

  // evaluate the kernel
  double wk; kernel_main(s, 1, 1, &wk, &wk, -1);

  // return the integrand
  return s * s * wk * result;
}


/**
 * tabulate values of the overlap integral.
 */
void dmb_init_overlap_table(void) {
  int i, j, index;
  double result, abserr, h_ratio, delta;
  double params[2];

  double t0 = my_second();
  PRINT_STATUS("precomputing the DM-b overlap integral table (%d ranks)...", NTask);

  // set up the integration workspace
  gsl_function F;
  gsl_integration_workspace *workspace;
  workspace = gsl_integration_workspace_alloc(GSLWORKSIZE);

  for (i = 0; i < DMB_OVERLAP_TABLE_LENGTH; ++i) {
    DMB_OverlapTable[i] = 0.0;
  }

  // distribute rows across MPI ranks
  for (i = ThisTask; i < DMB_OVERLAP_NUM_H_RATIO; i += NTask) {
    h_ratio = pow(10, -2.0 * ((double) i) / ((double) DMB_OVERLAP_NUM_H_RATIO - 1.0));

    for (j = 0; j < DMB_OVERLAP_NUM_DELTA; ++j) {
      index = i * DMB_OVERLAP_NUM_DELTA + j;

      delta = 2.0 * ((double) j) / ((double) DMB_OVERLAP_NUM_DELTA - 1.0);
      if (delta >= h_ratio + 1) {
        // no overlap, don't waste your time here or for any larger deltas
        break;
      }

      // register integrand and parameters
      params[0] = h_ratio;
      params[1] = delta;
      F.function = &overlap_s_integrand;
      F.params = params;

      // perform the integral
      gsl_integration_qag(&F, 0.0, +1.0, 0, 1e-6, GSLWORKSIZE, GSL_INTEG_GAUSS31, workspace, &result, &abserr);

      // store the result
      DMB_OverlapTable[index] = 2 * M_PI * result;

      if (result == 0.0) {
        // g_ij is monotonically decreasing as function of delta
        // so if we can't integrate it here, don't do larger deltas
        break;
      }
    }
  }

  gsl_integration_workspace_free(workspace);

  // combine partial results from all ranks
#ifndef DOUBLEPRECISION
  MPI_Allreduce(MPI_IN_PLACE, DMB_OverlapTable, DMB_OVERLAP_TABLE_LENGTH, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
#else
  MPI_Allreduce(MPI_IN_PLACE, DMB_OverlapTable, DMB_OVERLAP_TABLE_LENGTH, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  double timeall = timediff(t0, my_second());
  PRINT_STATUS("  ..finished! (%f sec)", timeall);
}


#endif // DM_DMB
