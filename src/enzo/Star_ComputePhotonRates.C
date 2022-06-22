/***********************************************************************
/
/  COMPUTE STELLAR PHOTON EMISSION RATES
/
/  written by: John Wise
/  date:       November, 2005
/  modified1:
/
/  ---------- SPECIES --------
/  0 : HI
/  1 : HeI
/  2 : HeII
/  3 : Lyman-Werner (H2)
/ 
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "phys_constants.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

float ReturnValuesFromSpectrumTable(float ColumnDensity, float dColumnDensity, int mode);

// KH 2021/6/28
// Need to be embedded into ENZO coding structure after testing.
// Define the spectrum function. L = Lo * (13.6/energy) ^ index  (eV/s/Hz)
// KH 2022/5/10
// L is multipled by exp(-energy/threshold)
double QSO_Luminosity_powlaw(const float &energy, const double &Lo, const float &index)
{
    double L_nu = 0.0;
    double nu_HI = 13.6;    // eV
    const double Temperature_threshold = 1.0e6; //K
    const double Boltzmann_constant = 8.617333262145e-5;   // eV / Kelvin
    // 1E6 K => 87eV (K_B * 1E6 ~ 86.17 eV)
    const double energy_threshold = Temperature_threshold * Boltzmann_constant; 
    L_nu = Lo * pow(energy / nu_HI, index);     //  eV/Hz/s
    //  eV/Hz/s
    // L_nu = Lo * pow(energy / nu_HI, index) * exp(-energy/energy_threshold);

    return L_nu;
}


// KH 2022/6/21 support blackbody radiation spectrum
// input : energy (eV);Lo (ev/(s*Hz))(power law);Lo_bbr (1/eV^2) (blackbody)
//         index (power index for power law, L ~ (E/E_HI)^(index)) 
//         frac_pow, frac_body : ratio between powlaw and bbr
double QSO_Luminosity_powlaw_bbr(const float &energy, const double &Lo, const double &Lo_bbr,
       const float &index, float frac_pow, float frac_body)
{
    double L_nu = 0.0;
    double L_nu_pow = 0.0;
    double L_nu_bbr = 0.0;
    // weighting between power law and black body radiation.
    // Lv = Lo * (frac_pow * "powlaw" + frac_body * "black body")
    // double frac_pow  = 0.9;  
    // double frac_body = 0.1;
    // Normalisation
    double temp = frac_pow + frac_body;
    frac_pow  = frac_pow  / temp;
    frac_body = frac_body / temp;
    double nu_HI = 13.6;    // eV
    const double Temperature_threshold = 1.0e6; //K
    const double Boltzmann_constant = 8.617333262145e-5;   // eV / Kelvin
    // 1E6 K => 87eV (K_B * 1E6 ~ 86.17 eV)
    const double energy_threshold = Temperature_threshold * Boltzmann_constant; 
    L_nu_pow = Lo * pow(energy / nu_HI, index);     //  eV/Hz/s
    //  eV/Hz/s
    // L_nu = Lo * pow(energy / nu_HI, index) * exp(-energy/energy_threshold);
    // KH 20/6/2022: add blackbody spectrum
    // eV/Hz/s
    L_nu_bbr = Lo_bbr * (pow(energy, 3.0)/expm1(energy/energy_threshold));
    L_nu = frac_pow * L_nu_pow + frac_body * L_nu_bbr;

    return L_nu;
}


int Star::ComputePhotonRates(const float TimeUnits, int &nbins, float E[], double Q[])
{


  int i;
  double L_UV, cgs_convert, _mass;
  float x, x2, EnergyFractionLW, MeanEnergy, XrayLuminosityFraction;
  float Mform, EnergyFractionHeI, EnergyFractionHeII;
  if (this->Mass < 0.1)  // Not "born" yet
    _mass = this->FinalMass;
  else
    _mass = this->Mass;
  x = log10((float)(_mass));
  x2 = x*x;

  switch(ABS(this->type)) {

    /* Luminosities from Schaerer (2002) */

  case PopIII:
    nbins = (PopIIIHeliumIonization &&
	     !RadiativeTransferHydrogenOnly) ? 3 : 1;
#ifdef TRANSFER    
    if (!RadiativeTransferOpticallyThinH2) nbins++;
#endif
    E[0] = 28.0;
    E[1] = 30.0;
    E[2] = 58.0;
    E[3] = 12.8;
    _mass = max(min((float)(_mass), 500), 5);
    if (_mass > 9 && _mass <= 500) {
      Q[0] = pow(10.0, 43.61 + 4.9*x   - 0.83*x2);
      Q[1] = pow(10.0, 42.51 + 5.69*x  - 1.01*x2);
      Q[2] = pow(10.0, 26.71 + 18.14*x - 3.58*x2);
      Q[3] = pow(10.0, 44.03 + 4.59*x  - 0.77*x2);
    } else if (_mass > 5 && _mass <= 9) {
      Q[0] = pow(10.0, 39.29 + 8.55*x);
      Q[1] = pow(10.0, 29.24 + 18.49*x);
      Q[2] = pow(10.0, 26.71 + 18.14*x - 3.58*x2);
      Q[3] = pow(10.0, 44.03 + 4.59*x  - 0.77*x2);
    } // ENDELSE
    else {
      for (i = 0; i < nbins; i++) Q[i] = 0.0;
    }
    break;

    /* Average energy from Schaerer (2003) */

  case PopII:
    nbins = (StarClusterHeliumIonization && 
	     !RadiativeTransferHydrogenOnly) ? 3 : 1;
#ifdef TRANSFER    
    if (!RadiativeTransferOpticallyThinH2 &&
	MultiSpecies > 1) nbins++;
#endif
    EnergyFractionLW   = 1.288;
    EnergyFractionHeI  = 0.2951;
    EnergyFractionHeII = 2.818e-4;
    E[0] = 21.62; // eV (good for a standard, low-Z IMF)
    E[1] = 30.0;
    E[2] = 60.0;
    E[3] = 12.8;
    Q[0] = StarClusterIonizingLuminosity * _mass;
    if (StarClusterHeliumIonization) {
      Q[1] = EnergyFractionHeI * Q[0];
      Q[2] = EnergyFractionHeII * Q[0];
      Q[0] *= 1.0 - EnergyFractionHeI - EnergyFractionHeII;
    } else {
      Q[1] = 0.0;
      Q[2] = 0.0;
    }
    Q[3] = EnergyFractionLW * Q[0];
    break;

    /* Approximation to the multi-color disk and power law of an
       accreting BH (Kuhlen & Madau 2004; Alvarez et al. 2009) */

  case BlackHole:
    nbins = 1;
    XrayLuminosityFraction = 0.43;
    EnergyFractionLW = 1.51e-3;
    MeanEnergy = 93.0;  // eV
    E[0] = 460.0;
    E[1] = 0.0;
    E[2] = 0.0;
    E[3] = 12.8;
    Q[0] = 1.12e66 * PopIIIBHLuminosityEfficiency * XrayLuminosityFraction *
      this->last_accretion_rate / E[0];
//    Below is wrong!
//    Q[0] = 3.54e58 * PopIIIBHLuminosityEfficiency * XrayLuminosityFraction *
//      this->DeltaMass / E[0];
    Q[1] = 0.0;
    Q[2] = 0.0;
    Q[3] = EnergyFractionLW * (E[0]/MeanEnergy) * Q[0];
    break;

    /* Average quasar SED by Sazonov et al.(2004), where associated 
       spectral temperature is 2 keV, for accreting massive BH */

  case MBH:
    {
        //KH adapting Gaussian Quadrature and QSO spectrum.
        const int number_bins = 20;    // MAX_ENERGY_BINS = 20
        nbins = number_bins;
        XrayLuminosityFraction = 1.0;
        double GQ_weight[number_bins] = {3.15087663e+14, 6.36525437e+14, 7.56561967e+14, 6.36525437e+14, 3.15087663e+14, 8.53601122e+14, 1.72440527e+15, 2.04959515e+15, 1.72440527e+15, 8.53601122e+14, 7.62203685e+15, 1.70856566e+16, 2.50465076e+16, 3.07832530e+16, 3.37850775e+16, 3.37850775e+16, 3.07832530e+16, 2.50465076e+16, 1.70856566e+16, 7.62203685e+15};  // 5(HI)+5(HeI)+10(HeII)
        /* Original ENZO setting
        E[0] = 2000.0; //2keV
        E[1] = 0.0;
        E[2] = 0.0;
        E[3] = 0.0;
        */
        E[0] = 14.11601085; 
        E[1] = 16.13841879;
        E[2] = 19.1;
        E[3] = 22.06158121;
        E[4] = 24.08398915;
        E[5] = 25.9979203;
        E[6] = 31.47680728;
        E[7] = 39.5;
        E[8] = 47.52319272; 
        E[9] = 53.0020797;
        E[10] = 66.73699332;
        E[11] = 118.19804023;
        E[12] = 205.97515611;
        E[13] = 322.29065766;
        E[14] = 456.81221253;
        E[15] = 597.58778747;
        E[16] = 732.10934234;
        E[17] = 848.42484389;
        E[18] = 936.20195977;
        E[19] = 987.66300668;
        // KH 2021: Use Eddington luminosity. (For long term, this should be achieved by adding a new module.)
        // Unit of Mass in star class is SolarMass
        // L_solar = 3.826 * 10 ^ 26 Joules/s = 2.388e45 eV/s
        double L_solar = 2.388e45;    // eV/s
        double L_edd   = 0.0;          // eV/s
        double L_total = 5.329276413387563e+58;  // eV/s Exact luminosity of QSO
        // double Lo_powlaw = 1.236708e43;    // eV/s/Hz for QSO
        double Lo_powlaw = 1.03140938373e+41;    // eV/s/Hz for 21cm
        double N_dot_photon = 0.0;                // Total emitted photon per second

        float  spectral_index = -1.73;
        // float  spectral_index = -0.5;     // for 21cm
        // float  spectral_index = 0.5;     // for 21cm
        float  reduced_factor = 0.01;   
        double Mass_Scaling_facter = 1e6;   // The Input Mass of MBHs has been reducded. Put the scaling back when calculate the luminosity. Default :1e8
        // L_edd = 3.2 * 10^4 * (Mass/Mass_solar) * L_solar
        L_edd = reduced_factor * 3.2 * 10000.0 * this->Mass * Mass_Scaling_facter * L_solar;   // eV/s
        // L_edd = 3.2 * 10000.0 * this->Mass * Mass_Scaling_facter * L_solar;   // eV/s
        double L_factor = Lo_powlaw * (L_edd / L_total);    // L_factor = Lo_powlaw * L_edd / L_total
        // KH 6/21/2022: for black body radiation
        const double Temperature_threshold = 1.0e6; //K
        const double Boltzmann_constant = 8.617333262145e-5;   // eV / Kelvin
        // 1E6 K => 87eV (K_B * 1E6 ~ 86.17 eV)
        const double energy_threshold = Temperature_threshold * Boltzmann_constant;
        double Lo_blackbody  = 0.0;    // 1/(eV^2) "bbr luminosity" at E = 13.6 
        double temp_integral = 0.0;

        N_dot_photon = 0.0;
        // Computing the total number of ionisation photon, based on power law 
        for(int j=0; j<nbins; j++)
        {
            N_dot_photon += (QSO_Luminosity_powlaw(E[j], Lo_powlaw, spectral_index) / E[j]) * GQ_weight[j]; // use a given luminosity
        }
        // Estimate Lo_blackbody (1/eV^2), expm1 = exp(x) - 1
        temp_integral = 0.0; 
        for(int j=0; j<nbins; j++)
        {
            temp_integral += (pow(E[j], 2.0) / expm1(E[j]/energy_threshold)) * GQ_weight[j];
        }
        Lo_blackbody = N_dot_photon / temp_integral; 
 
        for(int j=0; j<nbins; j++)
        {
            // Number of photon = (pow(13.6/energy, spectral_index) / energy) * (Lo_powlaw * L_edd/L_total) * GQ_weight    GQ_weight : weight computed by GQ.
            // Q[j] = (QSO_Luminosity_powlaw(E[j], L_factor, spectral_index) / E[j]) * GQ_weight[j]; // use eddington luminosity
            // Q[j] = (QSO_Luminosity_powlaw(E[j], Lo_powlaw, spectral_index) / E[j]) * GQ_weight[j]; // use a given luminosity
            // power law + blackboday radiation
            Q[j] = QSO_Luminosity_powlaw_bbr(E[j], Lo_powlaw, Lo_blackbody,
                    spectral_index, 0.9, 0.1);
        }
        for(int j=nbins; j<MAX_ENERGY_BINS; j++)
        {
            // Set the unused values to zero, following the original ENZO code.
            Q[j] = 0.0;
            E[j] = 0.0; 
        } 
        /* Original code
        // 1.99e33g/Ms * (3e10cm/s)^2 * 6.24e11eV/ergs = 1.12e66 eV/Ms (KH Ms=Msolar)
        // 1.12e66 eV/Msolar * accretion (Msolar/s) = eV/s 
        Q[0] = 1.12e66 * MBHFeedbackRadiativeEfficiency * XrayLuminosityFraction * this->last_accretion_rate / E[0];
        Q[0] = L_edd / E[0];
        Q[1] = 0.0;
        Q[2] = 0.0;
        Q[3] = 0.0;
        End */
        //KH testing
        fprintf(stderr, "function = %s, ", __FUNCTION__);
        fprintf(stderr, "ID = %lld, ", Identifier);
        fprintf(stderr, "star_type = %ld, ", type);
        fprintf(stderr, "L_edd = %g, ", L_edd);
        fprintf(stderr, "Q[0] = %g, ", Q[0]);
        fprintf(stderr, "Mass = %g (SolarMass). ", this->Mass);
        //fprintf(stderr, "MBHFeedbackRadiativeEfficiency = %g, ", MBHFeedbackRadiativeEfficiency);
        //fprintf(stderr, "XrayLuminosityFraction = %g, ", XrayLuminosityFraction);
        //fprintf(stderr, "this->last_accretion_rate = %g, ", this->last_accretion_rate);
        //fprintf(stderr, "last_accretion_rate = %g, ", last_accretion_rate);
        //fprintf(stderr, "last_accretion_rate = %p,  this->last_accretion_rate = %p, ", &last_accretion_rate, &this->last_accretion_rate);
        fprintf(stderr, "BirthTime = %"FSYM", LifeTime = %"FSYM".\n", BirthTime, LifeTime); 

#define NOT_HII_REGION_TEST
#ifdef HII_REGION_TEST
        Q[0] = 1.0e45 * MBHFeedbackRadiativeEfficiency * XrayLuminosityFraction / E[0];
#endif
    
//    fprintf(stdout, "star::ComputePhotonRates: this->last_accretion_rate = %g, Q[0]=%g\n", 
//    	    this->last_accretion_rate, Q[0]); 

#ifdef TRANSFER

        if (RadiativeTransferTraceSpectrum == TRUE) {
        nbins = 1;
        E[0] = ReturnValuesFromSpectrumTable(0.0, 0.0, 3); //##### mean energy if column density=0
        E[1] = 0.0;
        E[2] = 0.0;
        E[3] = 0.0;

        Q[0] = 1.12e66 * MBHFeedbackRadiativeEfficiency *
	    this->last_accretion_rate / E[0]; 
        Q[1] = 0.0;
        Q[2] = 0.0;
        Q[3] = 0.0;  

      //better check the initial mean energy when tracing spectrum
        if (MyProcessorNumber == ROOT_PROCESSOR)
	    fprintf(stdout, "star::CPP: check initial mean E of photon SED: E[0] = %g\n", E[0]); 
        }

#endif
    }
    break;

  case SimpleSource:
    nbins = 1;
    // radiating particle that ramps with time, independant of mass
    E[0] = 20.0;
    Q[0] = SimpleQ; // ramping done in StarParticleRadTransfer.C
    break;

  case NormalStar:
    nbins = 1;
    E[0] = 21.0;  // Good for [Z/H] > -1.3  (Schaerer 2003)
    // Calculate Delta(M_SF) for Cen & Ostriker star particles
#ifdef TRANSFER
    Mform = this->CalculateMassLoss(dtPhoton) / StarMassEjectionFraction;
    // units of Msun/(time in code units)
    L_UV = 4 * pi * StarEnergyToStellarUV * Mform * clight * clight / dtPhoton;
    cgs_convert = SolarMass / TimeUnits;
    Q[0] = cgs_convert * L_UV * eV_erg / E[0]; // ph/s
#else
    Q[0] = 0.0;
#endif
    break;

  default:
    ENZO_VFAIL("Star type = %"ISYM" not understood.\n", this->type)

  } // ENDSWITCH

  return SUCCESS;
}
