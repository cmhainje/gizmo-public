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

/*! Computes the momentum exchange rate per unit volume (which is filled into `out`).
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

/*! Computes the heat exchange rate per unit volume.
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

/*! Computes exchange rates and stores them in `out`. (out = [Pdot_x, Pdot_y, Pdot_z, Qdot]) */
void compute_exch_rates(int i, double out[4])
{
    int k;

    if (SphP[i].NumNgbDM == 0 || SphP[i].DM_Density <= 0) {
        for (k = 0; k < 4; k++) { out[k] = 0.0; }
        return;
    }

    // compute dV := V_DM - V_gas
    double dV[3];
    dV[0] = SphP[i].DM_Vx - P[i].Vel[0];
    dV[1] = SphP[i].DM_Vy - P[i].Vel[1];
    dV[2] = SphP[i].DM_Vz - P[i].Vel[2];

    // compute densities
    double rho_B = SphP[i].Density * All.cf_a3inv;
    double rho_DM = SphP[i].DM_Density * All.cf_a3inv;

    // compute temperatures (T_B computation follows galaxy_sf/sfr_eff.c:304-305)
    double u_B = SphP[i].InternalEnergyPred;
    double mu=1, ne=1, nh0=0, nHe0, nHepp, nhp, nHeII; // pull various known thermal properties, prepare to extract others //
    double T_B = ThermalProperties(u_B, rho_B, i, &mu, &ne, &nh0, &nhp, &nHe0, &nHeII, &nHepp); // get thermodynamic properties
    double T_DM = temperature_DM(SphP[i].DM_VelDisp);

    // compute 'microparticle' masses
    double m_B = mu * PROTONMASS_CGS; // TODO: gotta work out the units
    double m_DM = All.DMB_DarkMatterMass;

    // compute Pdot, Qdot
    double Pdot[3], Qdot;
    mom_exch_rate(dV, rho_DM, T_DM, m_DM, rho_B, T_B, m_B, Pdot);
    Qdot = heat_exch_rate(dV, rho_DM, T_DM, m_DM, rho_B, T_B, m_B);

    // fill output array
    for (k = 0; k < 3; k++) { out[k] = Pdot[k]; }
    out[3] = Qdot;
}



#endif
