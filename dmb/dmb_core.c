#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>

#include "../allvars.h"
#include "../proto.h"

/*! \file dmb_core.c
 *  \brief Fuctions and routines needed for the calculations of dark matter-baryon interactions
 *
 *  This file contains the functions and routines necesary for the computation of
 *  the momentum and heat exchange between dark matter and baryons due to interactions.
 *  Written by Connor Hainje, connor.hainje@nyu.edu, 2023.
 */

#ifdef DM_DMB

double HYPERG_ASYMP_FACTOR = 100;

/*! Computes the value of the momentum-exchange cross section at a given velocity.
 *  Units:
 *    velocity: physical, cgs (cm/s)
 *    return: physical, cgs (cm^2)
 */
double cross_section(double velocity)
{
    return All.DMB_InteractionCrossSection * pow(velocity, All.DMB_InteractionPowerScale);
}

/*! Implements the script A function (assuming a power-law cross section).
 *  Units:
 *    w: velocity, physical, CGS (cm/s)
 *    kT_over_m: physical (ergs/g) or (cm^2/s^2)
 */
double script_A(double w, double kT_over_m)
{
    // asymptotic limit, eq. 34
    if (w * w / kT_over_m > HYPERG_ASYMP_FACTOR) { return w * cross_section(w); }

    int n = All.DMB_InteractionPowerScale;
    double alpha = gsl_sf_hyperg_1F1(-0.5 * (n + 1), 2.5, -0.5 * w * w / kT_over_m);
    double sigma = All.DMB_InteractionCrossSection;
    double c = sqrt((1 << (5 + n)) / (9.0 * M_PI)) * gsl_sf_gamma(3.0 + 0.5 * n);
    //               ^^ (1 << k) computes 2^k as an int

    double out = c * sigma * pow(kT_over_m, 0.5 * (n + 1.0)) * alpha;

    if (isnan(out)) {
        printf("script_A returning NaN; inputs were w=%f, kT_over_m=%f\n", w, kT_over_m);
    }

    return out;
}

/*! Implements the script B function (assuming a power-law cross section).
 *  Units:
 *    w: velocity, physical, CGS (cm/s)
 *    kT_over_m: physical (ergs/g) or (cm^2/s^2)
 */
double script_B(double w, double kT_over_m)
{
    // asymptotic limit, eq. 34
    if (w * w / kT_over_m > HYPERG_ASYMP_FACTOR) { return w * w * w * cross_section(w); }

    int n = All.DMB_InteractionPowerScale;
    double beta = gsl_sf_hyperg_1F1(-0.5 * (n + 3), 1.5, -0.5 * w * w / kT_over_m);
    double sigma = All.DMB_InteractionCrossSection;
    double c = sqrt((1 << (5 + n)) / (9.0 * M_PI)) * gsl_sf_gamma(3.0 + 0.5 * n);
    //               ^^ (1 << k) computes 2^k as an int

    double out = 3.0 * c * sigma * pow(kT_over_m, 0.5 * (n + 3.0)) * beta;

    if (isnan(out)) {
        printf("script_B returning NaN; inputs were w=%f, kT_over_m=%f\n", w, kT_over_m);
    }

    return out;
}

/*! Computes the momentum exchange rate (B -> DM) per unit volume (which is filled into `out`).
 *  dV is the dark matter velocity minus the baryon velocity (in that order)
 *  rho_DM, kT_DM, m_DM are the mass density, temperature, and particle mass of the dark matter
 *  rho_B, kT_B, m_B are the same for baryonic matter
 * 
 *  Units:
 *    dV: velocity, physical, cgs (cm/s)
 *    rho_*: density, physical, cgs
 *    m_*: grams (physical)
 *    kT_*: ergs
 */
void mom_exch_rate(double dV[3], double rho_DM, double kT_DM, double m_DM, double rho_B, double kT_B, double m_B, double out[3])
{
    double dV_mag = sqrt(dV[0]*dV[0] + dV[1]*dV[1] + dV[2]*dV[2]);
    double v_th_2 = kT_B / m_B + kT_DM / m_DM;
    double A = script_A(dV_mag, v_th_2);
    double coeff = -(rho_DM * rho_B) / (m_DM + m_B) * A;

    int i;
    bool nan_detected = false;
    for (i = 0; i < 3; i++) { out[i] = coeff * dV[i]; nan_detected = nan_detected || isnan(out[i]); }
    if (nan_detected) {
        printf("mom_exch_rate returning NaN. inputs were:");
        printf("  rho_DM = %f\n", rho_DM);
        printf("  kT_DM = %f\n", kT_DM);
        printf("  m_DM = %f\n", m_DM);
        printf("  rho_B = %f\n", rho_B);
        printf("  kT_B = %f\n", kT_B);
        printf("  m_B = %f\n", m_B);
    }
}

/*! Computes the heat exchange rate (B -> DM) per unit volume.
 *  dV is the dark matter velocity minus the baryon velocity (in that order)
 *  rho_DM, kT_DM, m_DM are the mass density, temperature, and particle mass of the dark matter
 *  rho_B, kT_B, m_B are the same for baryonic matter
 * 
 *  Units:
 *    dV: velocity, physical, cgs (cm/s)
 *    rho_*: density, physical, cgs
 *    m_*: grams (physical)
 *    kT_*: ergs
 */
double heat_exch_rate(double dV[3], double rho_DM, double kT_DM, double m_DM, double rho_B, double kT_B, double m_B)
{
    double dV_mag = sqrt(dV[0]*dV[0] + dV[1]*dV[1] + dV[2]*dV[2]);
    double v_th_2 = kT_B / m_B + kT_DM / m_DM;
    double A = script_A(dV_mag, v_th_2);
    double B = script_B(dV_mag, v_th_2);
    double coeff = (rho_DM * rho_B) / (m_DM + m_B) / v_th_2;
    double out = coeff * (B * (kT_B - kT_DM) / (m_DM + m_B) + kT_DM / m_DM * A * dV_mag * dV_mag);

    if (isnan(out)) {
        printf("heat_exch_rate returning NaN. inputs were:\n");
        printf("  rho_DM = %e\n", rho_DM);
        printf("  kT_DM = %e\n", kT_DM);
        printf("  m_DM = %e\n", m_DM);
        printf("  rho_B = %e\n", rho_B);
        printf("  kT_B = %e\n", kT_B);
        printf("  m_B = %e\n", m_B);
    }

    return out;
}


/*! Computes the temperature of dark matter from its velocity dispersion.
 *  Assumes vel_disp is the 1D velocity dispersion and is given in code units.
 *  Returns temperature in ergs (e.g. returns kT).
 */
double temperature_DM(double vel_disp)
{
    double vd = vel_disp * UNIT_VEL_IN_CGS / All.cf_atime; // code -> phys
    return All.DMB_DarkMatterMass * vd * vd;
}

/*! Computes exchange rates and stores them in `accel` and `dUdt`. */
void compute_exch_rates_DM(int i, double accel[3], double *dUdt) {
    int k;

    // compute dV := v_DM (self) - v_gas (other) in [cgs]
    double dV[3]; for (k = 0; k < 3; k++) { dV[k] = (P[i].Vel[k] - P[i].DMB_V[k]) / All.cf_atime * UNIT_VEL_IN_CGS; }

    // densities
    double rho_DM = P[i].AGS_Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;
    double rho_gas = P[i].DMB_Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;

    // temperatures
    double kT_DM = P[i].DMB_MyTemp;
    double kT_gas = P[i].DMB_Temperature;

    // compute momentum, internal energy exchange rates per volume
    mom_exch_rate(dV, rho_DM, kT_DM, All.DMB_DarkMatterMass, rho_gas, kT_gas, P[i].DMB_GasMass, P[i].DMB_MomExch);
    P[i].DMB_HeatExch = heat_exch_rate(dV, rho_DM, kT_DM, All.DMB_DarkMatterMass, rho_gas, kT_gas, P[i].DMB_GasMass);

    // translate exchange rates into accel and d(spec energy)/dt in code units
    for (k = 0; k < 3; k++) { accel[k] = (P[i].DMB_MomExch[k] / rho_DM) / (UNIT_VEL_IN_CGS * UNIT_TIME_IN_CGS) * All.cf_atime; }
    *dUdt = (P[i].DMB_HeatExch / rho_DM) / (UNIT_SPECEGY_IN_CGS / UNIT_TIME_IN_CGS); // note: not sure if I need an All.cf_* factor here
}

/*! Computes exchange rates and stores them in `accel` and `dUdt`. */
void compute_exch_rates_gas(int i, double accel[3], double *dUdt) {
    int k;

    // compute dV := v_DM (other) - v_self (gas) in [cgs]
    double dV[3]; for (k = 0; k < 3; k++) { dV[k] = (P[i].DMB_V[k] - P[i].Vel[k]) / All.cf_atime * UNIT_VEL_IN_CGS; }

    // densities
    double rho_gas = SphP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;
    double rho_DM = P[i].DMB_Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;

    // temperatures
    double kT_gas = P[i].DMB_MyTemp;
    double kT_DM = P[i].DMB_Temperature;

    // compute B -> DM momentum, internal energy exchange rates per volume
    double Pdot_DM[3]; mom_exch_rate(dV, rho_DM, kT_DM, All.DMB_DarkMatterMass, rho_gas, kT_gas, P[i].DMB_MyMass, Pdot_DM);
    double Qdot_DM = heat_exch_rate(dV, rho_DM, kT_DM, All.DMB_DarkMatterMass, rho_gas, kT_gas, P[i].DMB_MyMass);

    // convert B -> DM into DM -> B
    for (k = 0; k < 3; k++) { P[i].DMB_MomExch[k] = -1 * Pdot_DM[k]; }
    P[i].DMB_HeatExch = (
        P[i].DMB_MomExch[0] * dV[0]
        + P[i].DMB_MomExch[1] * dV[1]
        + P[i].DMB_MomExch[2] * dV[2]
        - Qdot_DM
    );

    // translate exchange rates into accel and d(spec energy)/dt in code units
    for (k = 0; k < 3; k++) { accel[k] = (P[i].DMB_MomExch[k] / rho_gas) / (UNIT_VEL_IN_CGS * UNIT_TIME_IN_CGS) * All.cf_atime; }
    *dUdt = (P[i].DMB_HeatExch / rho_gas) / (UNIT_SPECEGY_IN_CGS / UNIT_TIME_IN_CGS); // note: not sure if I need an All.cf_* factor here
}

/*! Computes exchange rates and stores them in `accel` and `dUdt`. */
void compute_exch_rates(int i, double accel[3], double *dUdt) {
    if (P[i].DMB_NumNgb == 0) {
        int k; for (k = 0; k < 3; k++) { accel[k] = 0.; }
        *dUdt = 0.;
        return;
    }

    if (P[i].Type == 0) {
        compute_exch_rates_gas(i, accel, dUdt);
    } else if (P[i].Type == 1) {
        compute_exch_rates_DM(i, accel, dUdt);
    }

    // int k;

    // if (P[i].DMB_NumNgb == 0) {
    //     for (k = 0; k < 3; k++) { accel[k] = 0.; }
    //     *dUdt = 0.;
    //     return;
    // }

    // // compute dV := v_self - v_other in [cgs]
    // double dV[3]; for (k = 0; k < 3; k++) { dV[k] = (P[i].Vel[k] - P[i].DMB_V[k]) / All.cf_atime * UNIT_VEL_IN_CGS; }

    // // densities
    // double rho_self = ((P[i].Type == 0) ? SphP[i].Density : P[i].AGS_Density) * All.cf_a3inv * UNIT_DENSITY_IN_CGS;
    // double rho_other = P[i].DMB_Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;

    // // temperatures
    // double kT_self = P[i].DMB_MyTemp;
    // double kT_other = P[i].DMB_Temperature;

    // // compute momentum, internal energy exchange rates per volume
    // if (P[i].Type == 0) {
    //     mom_exch_rate(dV, rho_other, kT_other, All.DMB_DarkMatterMass, rho_self, kT_self, P[i].DMB_MyMass, P[i].DMB_MomExch);
    //     P[i].DMB_HeatExch = -1 * heat_exch_rate(dV, rho_other, kT_other, All.DMB_DarkMatterMass, rho_self, kT_self, P[i].DMB_MyMass);
    // } else {
    //     mom_exch_rate(dV, rho_self, kT_self, All.DMB_DarkMatterMass, rho_other, kT_other, P[i].DMB_GasMass, P[i].DMB_MomExch);
    //     P[i].DMB_HeatExch = heat_exch_rate(dV, rho_self, kT_self, All.DMB_DarkMatterMass, rho_other, kT_other, P[i].DMB_GasMass);
    // }

    // if (isnan(P[i].DMB_HeatExch)) {
    //     printf("  type = %d\n", P[i].Type);
    // }

    // // translate exchange rates into accel and d(spec energy)/dt in code units
    // for (k = 0; k < 3; k++) { accel[k] = (P[i].DMB_MomExch[k] / rho_self) / (UNIT_VEL_IN_CGS * UNIT_TIME_IN_CGS) * All.cf_atime; }
    // *dUdt = (P[i].DMB_HeatExch / rho_self) / (UNIT_SPECEGY_IN_CGS / UNIT_TIME_IN_CGS); // note: not sure if I need an All.cf_* factor here
}

/*! This function simply initializes some variables to prevent memory errors */
void dmb_init() {
    int i; for (i = 0; i < NumPart; i++) {
        P[i].DMB_Hsml = 0;
        P[i].DMB_NumNgb = 0;

        P[i].DMB_GasMass = 0;

        P[i].DMB_InternalEnergy = 0;
        P[i].DMB_HeatExch = 0;

        int k; for (k = 0; k < 3; k++) {
            P[i].DMB_MomExch[k] = 0;
        }
    }
}

#endif
