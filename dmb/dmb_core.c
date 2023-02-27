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
 *  Written by Connor Hainje, connor.hainje@nyu.edu, Dec 2022.
 */

#ifdef DM_DMB

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
 *    T_over_m: temp/mass, physical, CGS (K/g)
 */
double script_A(double w, double T_over_m)
{
    int n = All.DMB_InteractionPowerScale;
    double alpha = gsl_sf_hyperg_1F1(-0.5 * (n + 1), 2.5, -0.5 * w * w / T_over_m);
    double sigma = All.DMB_InteractionCrossSection;
    double c = sqrt((1 << (5 + n)) / (9.0 * M_PI)) * gsl_sf_gamma(3.0 + 0.5 * n);
    //               ^^ (1 << k) computes 2^k as an int

    double out = c * sigma * pow(T_over_m, 0.5 * (n + 1.0)) * alpha;

    if (isnan(out)) {
        printf("script_A returning NaN; inputs were w=%f, T_over_m=%f\n", w, T_over_m);
    }

    return out;
}

/*! Implements the script B function (assuming a power-law cross section).
 *  Units:
 *    w: velocity, physical, CGS (cm/s)
 *    T_over_m: temp/mass, physical, CGS (K/g)
 */
double script_B(double w, double T_over_m)
{
    int n = All.DMB_InteractionPowerScale;
    double beta = gsl_sf_hyperg_1F1(-0.5 * (n + 3), 1.5, -0.5 * w * w / T_over_m);
    double sigma = All.DMB_InteractionCrossSection;
    double c = sqrt((1 << (5 + n)) / (9.0 * M_PI)) * gsl_sf_gamma(3.0 + 0.5 * n);
    //               ^^ (1 << k) computes 2^k as an int

    double out = 3.0 * c * sigma * pow(T_over_m, 0.5 * (n + 3.0)) * beta;

    if (isnan(out)) {
        printf("script_B returning NaN; inputs were w=%f, T_over_m=%f\n", w, T_over_m);
    }

    return out;
}

/*! Computes the momentum exchange rate (B -> DM) per unit volume (which is filled into `out`).
 *  dV is the dark matter velocity minus the baryon velocity (in that order)
 *  rho_DM, T_DM, m_DM are the mass density, temperature, and particle mass of the dark matter
 *  rho_B, T_B, m_B are the same for baryonic matter
 * 
 *  Units:
 *    dV: velocity, physical, cgs (cm/s)
 *    rho_*: density, physical, cgs
 *    m_*: grams (physical)
 *    T_*: Kelvin
 */
void mom_exch_rate(double dV[3], double rho_DM, double T_DM, double m_DM, double rho_B, double T_B, double m_B, double out[3])
{
    double dV_mag = sqrt(dV[0]*dV[0] + dV[1]*dV[1] + dV[2]*dV[2]);
    double v_th_2 = T_B / m_B + T_DM / m_DM;
    double A = script_A(dV_mag, v_th_2);
    double coeff = -(rho_DM * rho_B) / (m_DM + m_B) * A;

    int i;
    bool nan_detected = false;
    for (i = 0; i < 3; i++) { out[i] = coeff * dV[i]; nan_detected = nan_detected || isnan(out[i]); }
    if (nan_detected) {
        printf("mom_exch_rate returning NaN. inputs were:");
        printf("  rho_DM = %f\n", rho_DM);
        printf("  T_DM = %f\n", T_DM);
        printf("  m_DM = %f\n", m_DM);
        printf("  rho_B = %f\n", rho_B);
        printf("  T_B = %f\n", T_B);
        printf("  m_B = %f\n", m_B);
    }
}

/*! Computes the heat exchange rate (B -> DM) per unit volume.
 *  dV is the dark matter velocity minus the baryon velocity (in that order)
 *  rho_DM, T_DM, m_DM are the mass density, temperature, and particle mass of the dark matter
 *  rho_B, T_B, m_B are the same for baryonic matter
 * 
 *  Units:
 *    dV: velocity, physical, cgs (cm/s)
 *    rho_*: density, physical, cgs
 *    m_*: grams (physical)
 *    T_*: Kelvin
 */
double heat_exch_rate(double dV[3], double rho_DM, double T_DM, double m_DM, double rho_B, double T_B, double m_B)
{
    double dV_mag = sqrt(dV[0]*dV[0] + dV[1]*dV[1] + dV[2]*dV[2]);
    double v_th_2 = T_B / m_B + T_DM / m_DM;
    double A = script_A(dV_mag, v_th_2);
    double B = script_B(dV_mag, v_th_2);
    double coeff = (rho_DM * rho_B) / (m_DM + m_B) / v_th_2;
    double out = coeff * (B * (T_B - T_DM) + T_DM / m_DM * A * dV_mag * dV_mag);

    if (isnan(out)) {
        printf("heat_exch_rate returning NaN. inputs were:");
        printf("  rho_DM = %f\n", rho_DM);
        printf("  T_DM = %f\n", T_DM);
        printf("  m_DM = %f\n", m_DM);
        printf("  rho_B = %f\n", rho_B);
        printf("  T_B = %f\n", T_B);
        printf("  m_B = %f\n", m_B);
    }

    return out;
}


/*! Computes the temperature of dark matter from its velocity dispersion.
 *  Assumes vel_disp is sigma^2 and is given in physical CGS units.
 *  Returns temperature in Kelvin.
 */
double temperature_DM(double vel_disp)
{
    return All.DMB_DarkMatterMass * vel_disp / 3 / BOLTZMANN_CGS;
}

/*! Computes exchange rates for a gas particle (DM -> B) and stores them in `out`. */
void compute_exch_rates_gas(int i, double pdot[3], double* qdot) {
    int k;

    if (P[i].DMB_NumNgbDM == 0 || P[i].DMB_DensityDM <= 0) {
        for (k = 0; k < 3; k++) { pdot[k] = 0.0; }
        *qdot = 0.;
        return;
    }

    // compute dV := V_DM - V_gas in [cm/s]
    double dV[3];
    for (k = 0; k < 3; k++) {
        dV[k] = (P[i].DMB_VDM[k] - P[i].Vel[k]) / All.cf_atime * UNIT_VEL_IN_CGS;
    }

    // compute densities in [cgs]
    double rho_B = SphP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;
    double rho_DM = P[i].DMB_DensityDM * All.cf_a3inv * UNIT_DENSITY_IN_CGS;

    // compute temperatures in [K]
    double u_B = SphP[i].InternalEnergyPred;
    double mu=1, ne=1, nh0=0, nHe0, nHepp, nhp, nHeII; // pull various known thermal properties, prepare to extract others //
    double T_B = ThermalProperties(u_B, rho_B, i, &mu, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // [K]

    double sigma_sq = P[i].DMB_VelDispDM * P[i].DMB_VelDispDM * UNIT_VEL_IN_CGS * UNIT_VEL_IN_CGS; // no co-moving factor; taken care of in hsml loop
    double T_DM = temperature_DM(sigma_sq); // [K]

    // compute 'microparticle' masses in [g]
    double m_B = mu * PROTONMASS_CGS;  // [g]
    double m_DM = All.DMB_DarkMatterMass;  // [g]

    // compute Pdot, Qdot
    double Pdot[3], Qdot;
    mom_exch_rate(dV, rho_DM, T_DM, m_DM, rho_B, T_B, m_B, Pdot);
    Qdot = heat_exch_rate(dV, rho_DM, T_DM, m_DM, rho_B, T_B, m_B);

    // fill output array
    for (k = 0; k < 3; k++) { pdot[k] = -1 * Pdot[k]; }
    *qdot = -1 * Qdot;
}

void compute_exch_rates_DM(int i, double pdot[3], double* qdot) {
    int k;

    if (P[i].DMB_NumNgbGas == 0 || P[i].DMB_DensityGas <= 0) {
        for (k = 0; k < 3; k++) { pdot[k] = 0.0; }
        *qdot = 0.;
        return;
    }

    // compute dV := V_DM - V_gas in [cgs]
    double dV[3];
    for (k = 0; k < 3; k++) {
        dV[k] = (P[i].Vel[k] - P[i].DMB_VGas[k]) / All.cf_atime * UNIT_VEL_IN_CGS;
    }

    // compute densities [cgs]
    double rho_B = P[i].DMB_DensityGas * All.cf_a3inv * UNIT_DENSITY_IN_CGS;
    double rho_DM = P[i].DMB_DensityDM * All.cf_a3inv * UNIT_DENSITY_IN_CGS;

    // compute temperatures [K]
    double T_B = P[i].DMB_TemperatureGas;
    double sigma_sq = P[i].DMB_VelDispDM * P[i].DMB_VelDispDM * UNIT_VEL_IN_CGS * UNIT_VEL_IN_CGS; // no comoving factor; taken care of in hsml loop
    double T_DM = temperature_DM(sigma_sq);

    // compute 'microparticle' masses [cgs]
    double m_B = P[i].DMB_MicroparticleMassGas * PROTONMASS_CGS;
    double m_DM = All.DMB_DarkMatterMass;

    // compute Pdot, Qdot
    double Pdot[3], Qdot;
    mom_exch_rate(dV, rho_DM, T_DM, m_DM, rho_B, T_B, m_B, Pdot);
    Qdot = heat_exch_rate(dV, rho_DM, T_DM, m_DM, rho_B, T_B, m_B);

    // fill output
    for (k = 0; k < 3; k++) { pdot[k] = Pdot[k]; }
    *qdot = Qdot;
}

/*! Computes exchange rates and stores them in `out`. (out = [Pdot_x, Pdot_y, Pdot_z, Qdot]) */
void compute_exch_rates(int i, double pdot[3], double* qdot)
{
    if (P[i].Type == 0) { compute_exch_rates_gas(i, pdot, qdot); }
    else if (P[i].Type == 1) { compute_exch_rates_DM(i, pdot, qdot); }
}

void compute_kicks_gas(int i, double v_kick[3], double *q_kick) {
    int k;
    double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);

    double v_coeff = dt / (P[i].DMB_DensityGas * All.cf_a3inv);
    double v_units = 1 / UNIT_VEL_IN_CGS / All.cf_atime; // phys -> code
    for (k = 0; k < 3; k++) { v_kick[k] = v_coeff * P[i].DMB_MomExch[k] * v_units; }

    double q_coeff = dt * P[i].Mass / (SphP[i].Density * All.cf_a3inv);
    double q_units = 1 / UNIT_ENERGY_IN_CGS; // phys -> code
    *q_kick = q_coeff * P[i].DMB_HeatExch * q_units;
}

void compute_kicks_DM(int i, double v_kick[3], double *q_kick) {
    int k;
    double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);

    double v_coeff = dt / (P[i].DMB_DensityDM * All.cf_a3inv);
    double v_units = 1 / UNIT_VEL_IN_CGS / All.cf_atime; // phys -> code
    for (k = 0; k < 3; k++) { v_kick[k] = v_coeff * P[i].DMB_MomExch[k] * v_units; }

    *q_kick = 0;  // not tracking DM internal energy at present
}

/*! Computes the velocity and internal energy kicks for particle i.
 *  The momentum and heat exchange rates are per unit time and unit volume
 *  The velocity kick is thus:
 *    dV := MomExchRate / M * dt * vol = MomExchRate * dt / density
 *  The energy kick is thus:
 *    dQ := HeatExchRate * dt * vol,    vol estimated by (M / density)
 *  The kicks returned are in code units.
 */
void compute_kicks(int i, double v_kick[3], double *q_kick) {
    if (P[i].Type == 0) { compute_kicks_gas(i, v_kick, q_kick); }
    else if (P[i].Type == 1) { compute_kicks_DM(i, v_kick, q_kick); }
}

#endif
