// Relevant physics constants
#ifndef __RUN1BANA_ANALYSIS_PHYSICS__
#define __RUN1BANA_ANALYSIS_PHYSICS__

//--------------------------------------------------------
// Livetime/muon stop normalization

double duty_cycle_1bb_ = 0.323                 ; // taken from https://github.com/Mu2e/Production/blob/main/JobConfig/ensemble/python/normalizations.py
double duty_cycle_2bb_ = 0.246                 ;

// SU2020 info
double livetime_su2020_        = 1.11e7;
double npot_su2020_            = 3.8e19;
double nmuons_per_pot_su2020_  = 0.00159;
double nmuons_su2020_          = npot_su2020_*nmuons_per_pot_su2020_;

// Run 1B defaults
double nmuons_per_pot_run1b_   = 9.561e-04;
double npot_per_s_1bb_         = duty_cycle_1bb_ * 1.6e7 / 1695.e-9; // duty factor * POT/cycle/1695.e-9 (1BB: 0.323 duty cycle, 1.6e7 POT / cycle)
double livetime_month_         = 2.6166667e6; // appox 1/12th of a year * duty factor: 3.14e7/12 (1BB)
double livetime_week_          = 6.0384616e5; // appox 1/12th of a year * duty factor: 3.14e7/12 (1BB)

// One month, 1BB
double npot_month_1bb_         = npot_per_s_1bb_*livetime_month_;
double nmuons_month_1bb_       = npot_month_1bb_*nmuons_per_pot_run1b_;

// One week, 1BB
double npot_week_1bb_          = npot_per_s_1bb_*livetime_week_;
double nmuons_week_1bb_        = npot_week_1bb_*nmuons_per_pot_run1b_;

// One week, 1% 1BB
double npot_week_1bb_1pc_      = 0.01*npot_per_s_1bb_*livetime_week_;
double nmuons_week_1bb_1pc_    = npot_week_1bb_1pc_*nmuons_per_pot_run1b_;

// Normalization used, defaulting to 1 week of running, 1% intensity
double livetime_       = livetime_week_        ;
double npot_           = npot_week_1bb_1pc_    ;
double nmuons_         = nmuons_week_1bb_1pc_  ;
double nmuons_per_pot_ = nmuons_per_pot_run1b_ ;
double signal_br_      = 1.e-9                 ; // example signal branching fraction for normalization

//--------------------------------------------------------
// Branching fractions

// MDS1g:
// pistops/POT 0.00223279
// pi time eff 0.06459918608102208
// pi surv prob 0.0008919234313257805
// BRRPC= 0.0215

double muon_capture_fraction_ = 0.609    ; // on aluminum
double muon_capture_carbon_   = 0.0701   ;
double ipa_nmuons_per_pot_    = 2.062e-08; // stopped muon / POT in IPA
double ipa_dio_frac_70_       = 2.538e-06; //fraction of carbon DIO spectrum above 70 MeV
// double pion_stop_rate_        = 0.00223  ; //Pion stops per POT;    SU2020: 0.00211; //N(pion) at ST per POT
// double pion_survive_frac_     = 0.064599 ; //Fraction above 350 ns; SU2020: 0.02423; //fraction of pions that survive 450 ns lifetime cut
double pion_stop_rate_        = 0.0002148; //Pion stops above 350 ns per POT;
double pion_survive_frac_     = 1.       ; //included in the number above
double rpc_frac_50_           = 0.9888   ; //fraction of RPC > 50 MeV, cutoff used in simulation
double rpc_br_                = 0.0215   ;
double rpc_int_br_            = 0.0069   ; //using BR(internal RPC) / BR(RPC) on hydrogen
double br_rmc_                = 1.40e-5  ; //For E > 57 MeV, relative to muon capture
double rmc_frac_57_           = 0.14800  ; //fraction of RMC > 57 (kmax = 90.1)
double rmc_frac_85_           = 0.0010641; //fraction of RMC > 85 (kmax = 90.1)
double rmc_int_br_            = 0.0069   ; //using BR(internal RPC) / BR(RPC) on hydrogen
double rmc_conv_              = 380711./1.e8; // fraction of RMC photons that convert
double pbar_stops_per_pot_    = 4.7e-18  ; //N(pbar at ST) / POT
double dio_frac_50_           = 8.766e-02; //DIO spectrum fraction above 50 MeV
double dio_frac_60_           = 2.735e-04; //DIO spectrum fraction above 60 MeV
double dio_frac_80_           = 5.671e-08; //DIO spectrum fraction above 80 MeV
double dio_frac_90_           = 7.262e-10; //DIO spectrum fraction above 90 MeV
double dio_frac_95_           = 3.637e-11; //DIO spectrum fraction above 95 MeV
double dio_frac_0_60_         = 0.999728 ; //DIO spectrum fraction  0 - 60
double dio_frac_60_80_        = 2.734e-04; //DIO spectrum fraction 60 - 80
double dio_frac_80_90_        = 5.950e-08; //DIO spectrum fraction 80 - 90

void init_physics(TString tag) {
  livetime_       = livetime_week_        ;
  npot_           = npot_week_1bb_1pc_    ;
  nmuons_         = nmuons_week_1bb_1pc_  ;
  nmuons_per_pot_ = nmuons_per_pot_run1b_ ;
  signal_br_      = 1.e-9                 ;
  if(tag.Contains("1bb_month")) {
    livetime_       = livetime_month_        ;
    npot_           = npot_month_1bb_    ;
    nmuons_         = nmuons_month_1bb_  ;
    nmuons_per_pot_ = nmuons_per_pot_run1b_ ;
  }
}

#endif
