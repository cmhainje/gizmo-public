#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>

#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/*! \file dmb_core.c
 *  \brief Fuctions and routines needed for the calculations of dark matter-baryon interactions
 *
 *  This file contains the functions and routines necesary for the computation of
 *  the momentum and heat exchange between dark matter and baryons due to interactions.
 *  Written by Connor Hainje, connor.hainje@nyu.edu, Dec 2022.
 */

#ifdef DM_DMB

/*! Computes the value of the momentum-exchange cross section at a given velocity. */
double cross_section(double velocity)
{
    return All.DMB_InteractionCrossSection * pow(velocity, All.DMB_InteractionPowerScale);
}

/*! Implements the script A function (assuming a power-law cross section). */
double script_A(double w, double T_over_m)
{
    int n = All.DMB_InteractionPowerScale;
    double alpha = gsl_sf_hyperg_1F1(-0.5 * (n + 1), 2.5, -0.5 * w * w / T_over_m);
    double sigma = cross_section(T_over_m);
    double c = sqrt((1 << (5 + n)) / (9.0 * M_PI)) * gsl_sf_gamma(3.0 + 0.5 * n);
    //               ^^ (1 << k) computes 2^k as an int

    return c * pow(sigma, 0.5 * (n + 1.0)) * alpha;
}

/*! Implements the script B function (assuming a power-law cross section). */
double script_B(double w, double T_over_m)
{
    int n = All.DMB_InteractionPowerScale;
    double beta = gsl_sf_hyperg_1F1(-0.5 * (n + 3), 1.5, -0.5 * w * w / T_over_m);
    double sigma = cross_section(T_over_m);
    double c = sqrt((1 << (5 + n)) / (9.0 * M_PI)) * gsl_sf_gamma(3.0 + 0.5 * n);
    //               ^^ (1 << k) computes 2^k as an int

    return 3.0 * c * pow(sigma, 0.5 * (n + 3.0)) * beta;
}

/*! Computes the momentum exchange rate (B -> DM) per unit volume (which is filled into `out`).
 *  dV is the dark matter velocity minus the baryon velocity (in that order)
 *  rho_DM, T_DM, m_DM are the mass density, temperature, and particle mass of the dark matter
 *  rho_B, T_B, m_B are the same for baryonic matter
 */
void mom_exch_rate(double dV[3], double rho_DM, double T_DM, double m_DM, double rho_B, double T_B, double m_B, double out[3])
{
    double dV_mag = sqrt(dV[0]*dV[0] + dV[1]*dV[1] + dV[2]*dV[2]);
    double v_th_2 = T_B / m_B + T_DM / m_DM;
    double A = script_A(dV_mag, v_th_2);
    double coeff = -(rho_DM * rho_B) / (m_DM + m_B) * A;

    int i;
    for (i = 0; i < 3; i++) { out[i] = coeff * dV[i]; }
}

/*! Computes the heat exchange rate (B -> DM) per unit volume.
 *  dV is the dark matter velocity minus the baryon velocity (in that order)
 *  rho_DM, T_DM, m_DM are the mass density, temperature, and particle mass of the dark matter
 *  rho_B, T_B, m_B are the same for baryonic matter
 */
double heat_exch_rate(double dV[3], double rho_DM, double T_DM, double m_DM, double rho_B, double T_B, double m_B)
{
    double dV_mag = sqrt(dV[0]*dV[0] + dV[1]*dV[1] + dV[2]*dV[2]);
    double v_th_2 = T_B / m_B + T_DM / m_DM;
    double A = script_A(dV_mag, v_th_2);
    double B = script_B(dV_mag, v_th_2);
    double coeff = (rho_DM * rho_B) / (m_DM + m_B) / v_th_2;
    return coeff * (B * (T_B - T_DM) + T_DM / m_DM * A * dV_mag * dV_mag);
}


/*! Computes the temperature of dark matter from its velocity dispersion.
 *  (Assumes vel_disp is sigma^2, not just sigma.)
 */
double temperature_DM(double vel_disp)
{
    // TODO: make sure units are correct (they probably aren't)
    return All.DMB_DarkMatterMass * vel_disp / 3;
}

/*! Computes exchange rates for a gas particle (DM -> B) and stores them in `out`. */
void compute_exch_rates_gas(int i, double pdot[3], double* qdot) {
    int k;

    if (P[i].DMB_NumNgbDM == 0 || P[i].DMB_DensityDM <= 0) {
        for (k = 0; k < 3; k++) { pdot[k] = 0.0; }
        *qdot = 0.;
        return;
    }

    // compute dV := V_DM - V_gas
    double dV[3];
    for (k = 0; k < 3; k++) { dV[k] = P[i].DMB_VDM[k] - P[i].Vel[k]; }

    // compute densities
    double rho_B = SphP[i].Density * All.cf_a3inv;
    double rho_DM = P[i].DMB_DensityDM * All.cf_a3inv;

    // compute temperatures
    double u_B = SphP[i].InternalEnergyPred;
    double mu=1, ne=1, nh0=0, nHe0, nHepp, nhp, nHeII; // pull various known thermal properties, prepare to extract others //
    double T_B = ThermalProperties(u_B, rho_B, i, &mu, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties
    double T_DM = temperature_DM(P[i].DMB_VelDispDM);

    // compute 'microparticle' masses
    double m_B = mu * PROTONMASS_CGS; // TODO: gotta work out the units
    double m_DM = All.DMB_DarkMatterMass;

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

    // compute dV := V_DM - V_gas
    double dV[3];
    for (k = 0; k < 3; k++) { dV[k] = P[i].Vel[k] - P[i].DMB_VGas[k]; }

    // compute densities
    double rho_B = P[i].DMB_DensityGas * All.cf_a3inv;
    double rho_DM = P[i].DMB_DensityDM * All.cf_a3inv;

    // compute temperatures
    double T_B = P[i].DMB_TemperatureGas;
    double T_DM = temperature_DM(P[i].DMB_VelDispDM);

    // compute 'microparticle' masses
    double m_B = P[i].DMB_MicroparticleMassGas * PROTONMASS_CGS; // TODO: gotta work out the units
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

#endif
