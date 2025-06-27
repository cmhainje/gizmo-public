#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_randist.h>

#include "../allvars.h"
#include "../proto.h"

/*! \file dmb_core.c
 *  \brief Fuctions and routines needed for the calculations of dark matter-baryon interactions
 *
 *  This file contains the functions and routines necesary for the computation of
 *  the momentum and heat exchange between dark matter and baryons due to interactions.
 *  Written by Connor Hainje, connor.hainje@nyu.edu, 2023-2025.
 */

#ifdef DM_DMB

double HYPERG_ASYMP_FACTOR = 100;

/*! Computes the value of the transfer cross section at a given velocity.
 *  Units:
 *    velocity: physical, cgs (cm/s)
 *    return: physical, cgs (cm^2)
 */
double transfer_cross_section(double velocity)
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
    if (w * w > HYPERG_ASYMP_FACTOR * kT_over_m)
        return w * transfer_cross_section(w);

    int n = All.DMB_InteractionPowerScale;
    double alpha = gsl_sf_hyperg_1F1(-0.5 * (n + 1), 2.5, -0.5 * w * w / kT_over_m);
    double sigma = All.DMB_InteractionCrossSection;
    double c = sqrt(pow(2., 5. + n) / (9.0 * M_PI)) * gsl_sf_gamma(3.0 + 0.5 * n);

    double out = c * sigma * pow(kT_over_m, 0.5 * (n + 1.0)) * alpha;

    if (isnan(out))
        printf("script_A returning NaN. inputs: w=%f, kT_over_m=%f\n", w, kT_over_m);

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
    if (w * w > HYPERG_ASYMP_FACTOR * kT_over_m)
        return w * w * w * transfer_cross_section(w);

    int n = All.DMB_InteractionPowerScale;
    double beta = gsl_sf_hyperg_1F1(-0.5 * (n + 3), 1.5, -0.5 * w * w / kT_over_m);
    double sigma = All.DMB_InteractionCrossSection;
    double c = sqrt(pow(2., 5. + n) / (9.0 * M_PI)) * gsl_sf_gamma(3.0 + 0.5 * n);

    double out = 3.0 * c * sigma * pow(kT_over_m, 0.5 * (n + 3.0)) * beta;

    if (isnan(out))
        printf("script_B returning NaN. inputs: w=%f, kT_over_m=%f\n", w, kT_over_m);

    return out;
}

/*! Computes the momentum exchange rate (B -> DM) per unit volume (which is filled into `out`).
 *  V_DMb is the dark matter velocity minus the baryon velocity (in that order)
 *  rho_DM, kT_DM, m_DM are the mass density, temperature, and particle mass of the dark matter
 *  rho_B, kT_B, m_B are the same for baryonic matter
 *
 *  Units:
 *    V_DMb: velocity, physical, cgs (cm/s)
 *    rho_*: density, physical, cgs
 *    m_*: grams (physical)
 *    kT_*: ergs
 */
void mom_exch_rate(double V_DMb[3], double rho_DM, double kT_DM, double m_DM, double rho_B, double kT_B, double m_B, double out[3])
{
    int k;

    if (rho_DM == 0.0)
    {
        for (k = 0; k < 3; k++)
            out[k] = 0.0;
        return;
    }

    double dV_mag = sqrt(V_DMb[0] * V_DMb[0] + V_DMb[1] * V_DMb[1] + V_DMb[2] * V_DMb[2]);
    double v_th_2 = kT_B / m_B + kT_DM / m_DM;
    double A = script_A(dV_mag, v_th_2);
    double coeff = -rho_B / (m_DM + m_B) * A;

    bool nan_detected = false;
    for (k = 0; k < 3; k++)
    {
        out[k] = coeff * V_DMb[k];
        nan_detected = nan_detected || isnan(out[k]);
    }
    if (nan_detected)
        printf(
            "mom_exch_rate returning NaN. inputs were:\n"
            "  rho_DM = %e\n"
            "  kT_DM = %e\n"
            "  m_DM = %e\n"
            "  rho_B = %e\n"
            "  kT_B = %e\n"
            "  m_B = %e\n",
            rho_DM, kT_DM, m_DM, rho_B, kT_B, m_B);
}

/*! Computes the heat exchange rate (B -> DM) per unit volume.
 *  V_DMb is the dark matter velocity minus the baryon velocity (in that order)
 *  rho_DM, kT_DM, m_DM are the mass density, temperature, and particle mass of the dark matter
 *  rho_B, kT_B, m_B are the same for baryonic matter
 *
 *  Units:
 *    V_DMb: velocity, physical, cgs (cm/s)
 *    rho_*: density, physical, cgs
 *    m_*: grams (physical)
 *    kT_*: ergs
 */
double heat_exch_rate(double V_DMb[3], double rho_DM, double kT_DM, double m_DM, double rho_B, double kT_B, double m_B)
{
    if (rho_DM == 0.0)
        return 0.0;

    double dV_mag = sqrt(V_DMb[0] * V_DMb[0] + V_DMb[1] * V_DMb[1] + V_DMb[2] * V_DMb[2]);
    double v_th_2 = kT_B / m_B + kT_DM / m_DM;
    if (v_th_2 == 0)
        return 0.0;

    double A = script_A(dV_mag, v_th_2);
    double B = script_B(dV_mag, v_th_2);
    double coeff = rho_B / (m_DM + m_B) / v_th_2;
    double out = coeff * (B * (kT_B - kT_DM) / (m_DM + m_B) + A * (kT_DM / m_DM) * dV_mag * dV_mag);

    if (isnan(out))
        printf(
            "heat_exch_rate returning NaN. inputs were:\n"
            "  rho_DM = %e\n"
            "  kT_DM = %e\n"
            "  m_DM = %e\n"
            "  rho_B = %e\n"
            "  kT_B = %e\n"
            "  m_B = %e\n",
            rho_DM, kT_DM, m_DM, rho_B, kT_B, m_B);

    return out;
}

/*! Computes the temperature of dark matter from its velocity dispersion.
 *  Assumes vel_disp is trace of the variance tensor and is given in code units of squared velocity.
 *  Returns temperature in ergs (e.g. returns kT).
 */
double temperature_DM(double vel_disp)
{
    double v_unit = UNIT_VEL_IN_CGS / All.cf_atime; // code -> phys velocity
    return All.DMB_DarkMatterMass * vel_disp * v_unit * v_unit / 3.0;
}

void print_everything(int i)
{
    // clang-format off
    printf(
        "  index %d, ID %d, type %d\n"
        "  xyz          = [%f, %f, %f]\n"
        "  vel          = [%f, %f, %f]\n"
        "  m_DM         = %e\n"
        "  AGS_Density  = %e\n"
        "  AGS_Hsml     = %e\n"
        "  AGS_NgbInt   = %d\n"
        "  AGS_NumNgb   = %e\n"
        "  AGS_VelMean  = [%f, %f, %f]\n"
        "  AGS_VelDisp  = %e\n"
        "  DMB_MyTemp   = %e\n"
        "  DMB_Hsml     = %e\n"
        "  DMB_NumNgb   = %e\n"
        "  DMB_V        = [%f, %f, %f]\n"
        "  DMB_Density  = %e\n"
        "  DMB_Temp     = %e\n"
        "  DMB_MyMass   = %e\n"
        "  DMB_GasMass  = %e\n"
        "  DMB_MomExch  = [%e, %e, %e]\n"
        "  DMB_HeatExch = %e\n"
        "  DMB_Accel    = [%e, %e, %e]\n"
        "  DMB_DtIntEgy = %e\n"
        "  DMB_MomExchd = %e\n"
        "  DMB_EgyExchd = %e\n",
        i, P[i].ID, P[i].Type,
        P[i].Pos[0], P[i].Pos[1], P[i].Pos[2],
        P[i].Vel[0], P[i].Vel[1], P[i].Vel[2],
        All.DMB_DarkMatterMass,
        P[i].AGS_Density,
        P[i].AGS_Hsml,
        P[i].AGS_NgbInt,
        P[i].AGS_NumNgb,
        P[i].AGS_VelMean[0], P[i].AGS_VelMean[1], P[i].AGS_VelMean[2],
        P[i].AGS_VelDisp,
        P[i].DMB_MyTemp,
        P[i].DMB_Hsml,
        P[i].DMB_NumNgb,
        P[i].DMB_V[0], P[i].DMB_V[1], P[i].DMB_V[2],
        P[i].DMB_Density,
        P[i].DMB_Temperature,
        P[i].DMB_MyMass,
        P[i].DMB_GasMass,
        P[i].DMB_MomExch[0], P[i].DMB_MomExch[1], P[i].DMB_MomExch[2],
        P[i].DMB_HeatExch,
        P[i].DMB_Accel[0], P[i].DMB_Accel[1], P[i].DMB_Accel[2],
        P[i].DMB_DtInternalEnergy,
        P[i].DMB_MomentumExchanged,
        P[i].DMB_EnergyExchanged);
    // clang-format on
}

/*! Computes exchange rates and stores them in `accel` and `dUdt`. */
void compute_exch_rates_DM(int i, double accel[3], double *dUdt)
{
    int k;

    // compute V_DMb := v_DM (self) - v_gas (other) in [cgs]
    double V_DMb[3];
    for (k = 0; k < 3; k++)
        V_DMb[k] = (P[i].AGS_VelMean[k] - P[i].DMB_V[k]) / All.cf_atime * UNIT_VEL_IN_CGS;

    // densities
    double rho_DM = P[i].AGS_Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;
    double rho_gas = P[i].DMB_Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;

    if (rho_DM == 0.0 || rho_gas == 0.0)
    {
        for (k = 0; k < 3; k++)
            accel[k] = 0.0;
        *dUdt = 0.0;
        return;
    }

    // temperatures (already in physical units)
    double kT_DM = P[i].DMB_MyTemp;
    double kT_gas = P[i].DMB_Temperature;

    // compute momentum, internal energy exchange rates per volume
    mom_exch_rate(V_DMb, rho_DM, kT_DM, All.DMB_DarkMatterMass, rho_gas, kT_gas, P[i].DMB_GasMass, P[i].DMB_MomExch);
    P[i].DMB_HeatExch = heat_exch_rate(V_DMb, rho_DM, kT_DM, All.DMB_DarkMatterMass, rho_gas, kT_gas, P[i].DMB_GasMass);

    // translate exchange rates into accel and d(spec energy)/dt in code units
    for (k = 0; k < 3; k++)
        accel[k] = P[i].DMB_MomExch[k] / (UNIT_VEL_IN_CGS / UNIT_TIME_IN_CGS) * All.cf_atime;
    *dUdt = P[i].DMB_HeatExch / (UNIT_SPECEGY_IN_CGS / UNIT_TIME_IN_CGS);
    // TODO: do I need an All.cf_* factor here?

#ifdef DM_DMB_NO_MOM
    for (k = 0; k < 3; k++)
        accel[k] = 0.0;
#endif
#ifdef DM_DMB_NO_HEAT
    *dUdt = 0.0;
#endif

    double accel_mag = sqrt(accel[0] * accel[0] + accel[1] * accel[1] + accel[2] * accel[2]);
    if ((accel_mag > 1e10) || (abs(*dUdt) > 1e10))
    {
        printf("Crazy behavior on ID %d (type %d):\n", i, P[i].Type);
        print_everything(i);
    }

    if ((GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i) > 0) && (GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i) < 1e-9))
    {
        printf("Tiny timestep for ID %d (type %d):\n", i, P[i].Type);
        print_everything(i);
    }

    // check for NaNs
    bool nan_detected = isnan(*dUdt);
    for (k = 0; k < 3; k++)
        nan_detected = nan_detected || isnan(accel[k]);

    if (nan_detected)
    {
        printf("compute_exch_rates_DM returning NaN on ID %d (type %d):\n", i, P[i].Type);
        print_everything(i);
    }
}

/*! Computes exchange rates and stores them in `accel` and `dUdt`. */
void compute_exch_rates_gas(int i, double accel[3], double *dUdt)
{
    int k;

    // compute V_DMb := v_DM (other) - v_self (gas) in [cgs]
    double V_DMb[3];
    for (k = 0; k < 3; k++)
        V_DMb[k] = (P[i].DMB_V[k] - P[i].Vel[k]) / All.cf_atime * UNIT_VEL_IN_CGS;

    // densities [g cm^-3]
    double rho_gas = SphP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;
    double rho_DM = P[i].DMB_Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;

    if (rho_DM == 0.0 || rho_gas == 0.0)
    {
        for (k = 0; k < 3; k++)
            accel[k] = 0.0;
        *dUdt = 0.0;
        return;
    }

    // temperatures [erg]
    double kT_gas = P[i].DMB_MyTemp;
    double kT_DM = P[i].DMB_Temperature;

    // compute B -> DM momentum, internal energy exchange rates per volume
    double Vdot_DM[3];
    mom_exch_rate(V_DMb, rho_DM, kT_DM, All.DMB_DarkMatterMass, rho_gas, kT_gas, P[i].DMB_MyMass, Vdot_DM);
    double Udot_DM = heat_exch_rate(V_DMb, rho_DM, kT_DM, All.DMB_DarkMatterMass, rho_gas, kT_gas, P[i].DMB_MyMass);

    // convert B -> DM into DM -> B
    double rho_ratio = rho_DM / rho_gas;
    for (k = 0; k < 3; k++)
        P[i].DMB_MomExch[k] = -rho_ratio * Vdot_DM[k];

    P[i].DMB_HeatExch = (P[i].DMB_MomExch[0] * V_DMb[0] + P[i].DMB_MomExch[1] * V_DMb[1] + P[i].DMB_MomExch[2] * V_DMb[2] - rho_ratio * Udot_DM);

    // translate exchange rates into accel and d(spec energy)/dt in code units
    for (k = 0; k < 3; k++)
        accel[k] = P[i].DMB_MomExch[k] / (UNIT_VEL_IN_CGS / UNIT_TIME_IN_CGS) * All.cf_atime;
    *dUdt = P[i].DMB_HeatExch / (UNIT_SPECEGY_IN_CGS / UNIT_TIME_IN_CGS);
    // TODO: do I need an All.cf_* factor here?

#ifdef DM_DMB_NO_MOM
    for (k = 0; k < 3; k++)
        accel[k] = 0.0;
#endif
#ifdef DM_DMB_NO_HEAT
    *dUdt = 0.0;
#endif

    double accel_mag = sqrt(accel[0] * accel[0] + accel[1] * accel[1] + accel[2] * accel[2]);
    if ((accel_mag > 1e10) || (abs(*dUdt) > 1e10))
    {
        printf("Crazy acceleration on ID %d (type %d):\n", i, P[i].Type);
        print_everything(i);
    }

    if ((GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i) > 0) && (GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i) < 1e-9))
    {
        printf("Tiny timestep for ID %d (type %d):\n", i, P[i].Type);
        print_everything(i);
    }

    // check for NaNs
    bool nan_detected = isnan(*dUdt);
    for (k = 0; k < 3; k++)
        nan_detected = nan_detected || isnan(accel[k]);

    if (nan_detected)
    {
        printf("compute_exch_rates_gas returning NaN on ID %d (type %d):\n", i, P[i].Type);
        print_everything(i);
    }
}

/*! Computes exchange rates and stores them in `DMB_Accel` and `DMB_DtInternalEnergy`. */
void compute_exch_rates(int i)
{
    if (P[i].DMB_NgbInt == 0)
    {
        int k;
        for (k = 0; k < 3; k++)
            P[i].DMB_Accel[k] = 0.;
        P[i].DMB_DtInternalEnergy = 0.;
        return;
    }

    if (P[i].Type == 0)
        compute_exch_rates_gas(i, P[i].DMB_Accel, &P[i].DMB_DtInternalEnergy);
    else if (P[i].Type == 1)
        compute_exch_rates_DM(i, P[i].DMB_Accel, &P[i].DMB_DtInternalEnergy);
}

/*! This function simply initializes some variables to prevent memory errors */
void dmb_init()
{
    int i;
    for (i = 0; i < NumPart; i++)
    {
        P[i].DMB_Hsml = 0;
        P[i].DMB_NumNgb = 0;
        P[i].DMB_NgbInt = 0;

        P[i].DMB_GasMass = 0;

        P[i].DMB_InternalEnergy = 0;
        // P[i].DMB_LastEnergyExchanged = 0;
        // P[i].DMB_EnergyError = 0;
        P[i].DMB_HeatExch = 0;

        int k;
        for (k = 0; k < 3; k++)
        {
            P[i].DMB_MomExch[k] = 0;
        }

        P[i].DMB_MomentumExchanged = 0;
        P[i].DMB_EnergyExchanged = 0;
    }
}

#endif
