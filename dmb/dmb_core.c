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
    if (w * w > HYPERG_ASYMP_FACTOR * kT_over_m) { return w * cross_section(w); }

    int n = All.DMB_InteractionPowerScale;
    double alpha = gsl_sf_hyperg_1F1(-0.5 * (n + 1), 2.5, -0.5 * w * w / kT_over_m);
    double sigma = All.DMB_InteractionCrossSection;
    double c = sqrt(pow(2., 5. + n) / (9.0 * M_PI)) * gsl_sf_gamma(3.0 + 0.5 * n);

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
    if (w * w > HYPERG_ASYMP_FACTOR * kT_over_m) { return w * w * w * cross_section(w); }

    int n = All.DMB_InteractionPowerScale;
    double beta = gsl_sf_hyperg_1F1(-0.5 * (n + 3), 1.5, -0.5 * w * w / kT_over_m);
    double sigma = All.DMB_InteractionCrossSection;
    double c = sqrt(pow(2., 5. + n) / (9.0 * M_PI)) * gsl_sf_gamma(3.0 + 0.5 * n);

    double out = 3.0 * c * sigma * pow(kT_over_m, 0.5 * (n + 3.0)) * beta;

    if (isnan(out)) {
        printf("script_B returning NaN; inputs were w=%f, kT_over_m=%f\n", w, kT_over_m);
    }

    return out;
}

/*! Computes the momentum exchange rate (2 -> 1) per unit volume (which is filled into `out`).
 *  dV is the velocity of species 1 minus species 2
 *  rho_1, kT_1, m_1 are the mass density, temperature, and particle mass of species 1
 *  rho_2, kT_2, m_2 are the same for species 2
 *
 *  Units:
 *    dV: velocity, physical, cgs (cm/s)
 *    rho_*: density, physical, cgs
 *    kT_*: ergs
 *    m_*: grams (physical)
 */
void mom_exch_rate(double dV[3], double rho_1, double kT_1, double m_1, double rho_2, double kT_2, double m_2, double out[3])
{
    double dV_mag = sqrt(dV[0]*dV[0] + dV[1]*dV[1] + dV[2]*dV[2]);
    double v_th_2 = kT_2 / m_2 + kT_1 / m_1;
    double A = script_A(dV_mag, v_th_2);
    double coeff = -(rho_1 * rho_2) / (m_1 + m_2) * A;

    int i;
    bool nan_detected = false;
    for (i = 0; i < 3; i++) { out[i] = coeff * dV[i]; nan_detected = nan_detected || isnan(out[i]); }
    if (nan_detected) {
        printf("mom_exch_rate returning NaN. inputs were:");
        printf("  rho_1 = %e\n", rho_1);
        printf("  kT_1  = %e\n", kT_1);
        printf("  m_1   = %e\n", m_1);
        printf("  rho_2 = %e\n", rho_2);
        printf("  kT_2  = %e\n", kT_2);
        printf("  m_2   = %e\n", m_2);
    }
}

/*! Computes the heat exchange rate (2 -> 1) per unit volume.
 *  dV is V_1 - V_2: mean velocity of species 1 minus 2
 *  rho_1, kT_1, m_1 are the mass density, temperature, and particle mass of species 1
 *  rho_2, kT_2, m_2 are the same for species 2
 *  correction is 3 (kT_DM / m_DM) / N_DM_neighbors
 *
 *  Units:
 *    dV: velocity, physical, cgs (cm/s)
 *    rho_*: density, physical, cgs
 *    kT_*: ergs
 *    m_*: grams (physical)
 *    correction: squared velocity, physical, cgs
 */
double heat_exch_rate(double dV[3], double rho_1, double kT_1, double m_1, double rho_2, double kT_2, double m_2, double correction)
{
    double dV_mag = sqrt(dV[0]*dV[0] + dV[1]*dV[1] + dV[2]*dV[2]);
    double v_th_2 = kT_2 / m_2 + kT_1 / m_1;
    double A = script_A(dV_mag, v_th_2);
    double B = script_B(dV_mag, v_th_2);
    double coeff = (rho_1 * rho_2) / (m_1 + m_2) / v_th_2;
    double out = coeff * (B * (kT_2 - kT_1) / (m_1 + m_2) + kT_1 / m_1 * A * fmax(dV_mag * dV_mag - correction, 0.0));
    if (v_th_2 == 0) out = 0;

    if (isnan(out)) {
        printf("heat_exch_rate returning NaN. inputs were:\n");
        printf("  rho_1 = %e\n", rho_1);
        printf("  kT_1  = %e\n", kT_1);
        printf("  m_1   = %e\n", m_1);
        printf("  rho_2 = %e\n", rho_2);
        printf("  kT_2  = %e\n", kT_2);
        printf("  m_2   = %e\n", m_2);
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
    return All.DMB_DarkMatterMass * vd * vd / 3.0;
}


void print_everything(int i) {
    if (P[i].Type == 0)
        printf(
            "  index %d, ID %d, type %d\n"
            "  xyz          = [%f, %f, %f]\n"
            "  vel          = [%f, %f, %f]\n"
            "  m_DM         = %e\n"
            "  SphP.Density = %e\n"
            "  DMB_MyTemp   = %e\n"
            "  DMB_Hsml     = %e\n"
            "  DMB_NumNgb   = %f\n"
            "  DMB_V        = [%f, %f, %f]\n"
            "  DMB_Density  = %e\n"
            "  DMB_Temp     = %e\n"
            "  DMB_MyMass   = %e\n"
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
            SphP[i].Density,
            P[i].DMB_MyTemp,
            P[i].DMB_Hsml,
            P[i].DMB_NumNgb,
            P[i].DMB_V[0], P[i].DMB_V[1], P[i].DMB_V[2],
            P[i].DMB_Density,
            P[i].DMB_Temperature,
            P[i].DMB_MyMass,
            P[i].DMB_MomExch[0], P[i].DMB_MomExch[1], P[i].DMB_MomExch[2],
            P[i].DMB_HeatExch,
            P[i].DMB_Accel[0], P[i].DMB_Accel[1], P[i].DMB_Accel[2],
            P[i].DMB_DtInternalEnergy,
            P[i].DMB_MomentumExchanged,
            P[i].DMB_EnergyExchanged);
    else if (P[i].Type == 1)
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
            P[i].DMB_GasMass,
            P[i].DMB_MomExch[0], P[i].DMB_MomExch[1], P[i].DMB_MomExch[2],
            P[i].DMB_HeatExch,
            P[i].DMB_Accel[0], P[i].DMB_Accel[1], P[i].DMB_Accel[2],
            P[i].DMB_DtInternalEnergy,
            P[i].DMB_MomentumExchanged,
            P[i].DMB_EnergyExchanged);
}



/*! Computes exchange rates and stores them in `accel` and `dUdt`. */
void compute_exch_rates_DM(int i, double accel[3], double *dUdt) {
    int k;

    // OLD: compute dV := v_DM (self) - v_gas (other) in [cgs]
    // compute dV := v_self - v_other in [cgs]
    double dV[3]; for (k = 0; k < 3; k++) { dV[k] = (P[i].AGS_VelMean[k] - P[i].DMB_V[k]) / All.cf_atime * UNIT_VEL_IN_CGS; }

    // densities
    double rho_DM = P[i].AGS_Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;
    double rho_gas = P[i].DMB_Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;

    if (rho_DM == 0.0 || rho_gas == 0.0) {
        for (k = 0; k < 3; k++) { accel[k] = 0.0; }
        *dUdt = 0.0;
        return;
    }

    // temperatures (already in physical units)
    double kT_DM = P[i].DMB_MyTemp;
    double kT_gas = P[i].DMB_Temperature;

    // compute momentum, internal energy exchange rates per volume
    mom_exch_rate(dV, rho_DM, kT_DM, All.DMB_DarkMatterMass, rho_gas, kT_gas, P[i].DMB_GasMass, P[i].DMB_MomExch);

    // double corr = 3. * kT_DM / All.DMB_DarkMatterMass / P[i].AGS_NgbInt;
    double corr = 0.0;
    P[i].DMB_HeatExch = heat_exch_rate(dV, rho_DM, kT_DM, All.DMB_DarkMatterMass, rho_gas, kT_gas, P[i].DMB_GasMass, corr);

    // translate exchange rates into accel and d(spec energy)/dt in code units
    for (k = 0; k < 3; k++) { accel[k] = (P[i].DMB_MomExch[k] / rho_DM) / (UNIT_VEL_IN_CGS / UNIT_TIME_IN_CGS) * All.cf_atime; }
    *dUdt = (P[i].DMB_HeatExch / rho_DM) / (UNIT_SPECEGY_IN_CGS / UNIT_TIME_IN_CGS); // note: not sure if I need an All.cf_* factor here

#ifdef DM_DMB_NO_MOM
    for (k = 0; k < 3; k++) { accel[k] = 0.0; }
#endif
#ifdef DM_DMB_NO_HEAT
    *dUdt = 0.0;
#endif

    double accel_mag = sqrt(accel[0]*accel[0] + accel[1]*accel[1] + accel[2]*accel[2]);
    if ((accel_mag > 1e10) || (abs(*dUdt) > 1e10)) {
        printf("Crazy behavior on ID %d (type %d):\n", P[i].ID, P[i].Type);
        print_everything(i);

        savepositions(999);
        endrun(1234);
    }

    if ((GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i) > 0) && (GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i) < 1e-9)) {
        printf("Tiny timestep for ID %d (type %d):\n", P[i].ID, P[i].Type);
        print_everything(i);

        savepositions(999);
        endrun(1234);
    }


    // check for NaNs
    bool nan_detected = isnan(*dUdt);
    for (k = 0; k < 3; k++) { nan_detected = nan_detected || isnan(accel[k]); }
    if (nan_detected) {
        printf("compute_exch_rates_DM returning NaN. inputs were:\n");
        print_everything(i);
    }
}

/*! Computes exchange rates and stores them in `accel` and `dUdt`. */
void compute_exch_rates_gas(int i, double accel[3], double *dUdt) {
    int k;

    // OLD: compute dV := v_DM (other) - v_self (gas) in [cgs]
    // compute dV := v_self - v_other in [cgs]
    double dV[3]; for (k = 0; k < 3; k++) { dV[k] = (P[i].Vel[k] - P[i].DMB_V[k]) / All.cf_atime * UNIT_VEL_IN_CGS; }

    // densities [g cm^-3]
    double rho_gas = SphP[i].Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;
    double rho_DM = P[i].DMB_Density * All.cf_a3inv * UNIT_DENSITY_IN_CGS;

    if (rho_DM == 0.0 || rho_gas == 0.0) {
        for (k = 0; k < 3; k++) { accel[k] = 0.0; }
        *dUdt = 0.0;
        return;
    }

    // temperatures [erg]
    double kT_gas = P[i].DMB_MyTemp;
    double kT_DM = P[i].DMB_Temperature;

    // compute energy exchange rates per volume
    mom_exch_rate(dV, rho_gas, kT_gas, P[i].DMB_MyMass, rho_DM, kT_DM, All.DMB_DarkMatterMass, P[i].DMB_MomExch);

    // double corr = 3. * kT_DM / All.DMB_DarkMatterMass / P[i].DMB_NgbInt;
    double corr = 0.0;
    P[i].DMB_HeatExch = heat_exch_rate(dV, rho_gas, kT_gas, P[i].DMB_MyMass, rho_DM, kT_DM, All.DMB_DarkMatterMass, corr);


    // translate exchange rates into accel and d(spec energy)/dt in code units
    for (k = 0; k < 3; k++) { accel[k] = (P[i].DMB_MomExch[k] / rho_gas) / (UNIT_VEL_IN_CGS / UNIT_TIME_IN_CGS) * All.cf_atime; }
    *dUdt = (P[i].DMB_HeatExch / rho_gas) / (UNIT_SPECEGY_IN_CGS / UNIT_TIME_IN_CGS); // note: not sure if I need an All.cf_* factor here

#ifdef DM_DMB_NO_MOM
    for (k = 0; k < 3; k++) { accel[k] = 0.0; }
#endif
#ifdef DM_DMB_NO_HEAT
    *dUdt = 0.0;
#endif

    double accel_mag = sqrt(accel[0]*accel[0] + accel[1]*accel[1] + accel[2]*accel[2]);
    if ((accel_mag > 1e10) || (abs(*dUdt) > 1e10)) {
        printf("Crazy behavior on ID %d (type %d):\n", P[i].ID, P[i].Type);
        print_everything(i);

        savepositions(999);
        endrun(1234);
    }

    if ((GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i) > 0) && (GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i) < 1e-9)) {
        printf("Tiny timestep for ID %d (type %d):\n", P[i].ID, P[i].Type);
        print_everything(i);

        savepositions(999);
        endrun(1234);
    }

    // check for NaNs
    bool nan_detected = isnan(*dUdt);
    for (k = 0; k < 3; k++) { nan_detected = nan_detected || isnan(accel[k]); }
    if (nan_detected) {
        printf("compute_exch_rates_gas returning NaN. inputs were:\n");
        print_everything(i);
    }
}

/*! Computes exchange rates and stores them in `accel` and `dUdt`. */
void compute_exch_rates(int i) {
    if (P[i].DMB_NgbInt == 0) {
        int k; for (k = 0; k < 3; k++) { P[i].DMB_Accel[k] = 0.; }
        P[i].DMB_DtInternalEnergy = 0.;
        return;
    }

    if (P[i].Type == 0)
        compute_exch_rates_gas(i, P[i].DMB_Accel, &P[i].DMB_DtInternalEnergy);
    else if (P[i].Type == 1)
        compute_exch_rates_DM(i, P[i].DMB_Accel, &P[i].DMB_DtInternalEnergy);
}

/*! This function simply initializes some variables to prevent memory errors */
void dmb_init() {
    int i; for (i = 0; i < NumPart; i++) {
        P[i].DMB_Hsml = 0;
        P[i].DMB_NumNgb = 0;
        P[i].DMB_NgbInt = 0;

        P[i].DMB_GasMass = 0;

        P[i].DMB_InternalEnergy = 0;
        // P[i].DMB_LastEnergyExchanged = 0;
        // P[i].DMB_EnergyError = 0;
        P[i].DMB_HeatExch = 0;

        int k; for (k = 0; k < 3; k++) {
            P[i].DMB_MomExch[k] = 0;
        }

        P[i].DMB_MomentumExchanged = 0;
        P[i].DMB_EnergyExchanged = 0;
    }
}

#endif
