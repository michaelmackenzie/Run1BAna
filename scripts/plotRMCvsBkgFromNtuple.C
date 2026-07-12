// Plot RMC vs. Bkg

#include "Run1BAna/scripts/plotSigvsBkgFromNtuple.C"

//------------------------------------------------------------------------------
void plotRMCvsBkgFromNtuple(const char* tag = "v40") {

  auto datasets = getDatasets(tag);
  if(datasets.empty()) {
    Error(__func__, "No datasets found for tag %s", tag);
    return;
  }

  TString sig_file = datasets.at("fgam").fileName();
  TString bkg_file = datasets.at("mnbs").fileName();
  TString csm_file = datasets.at("csms").fileName();
  TString neu_file = datasets.at("neut").fileName();
  TFile* f_sig = TFile::Open(sig_file, "READ");
  TFile* f_bkg = TFile::Open(bkg_file, "READ");
  TFile* f_csm = TFile::Open(csm_file, "READ");
  TFile* f_neu = TFile::Open(neu_file, "READ");
  if (!f_sig || f_sig->IsZombie() ||
      !f_bkg || f_bkg->IsZombie() ||
      !f_csm || f_csm->IsZombie() ||
      !f_neu || f_neu->IsZombie()
      ) {
    Error(__func__, "Could not open files!");
    return;
  }

  // General info
  const double onspill_time   = livetime_week_*duty_cycle_1bb_;
  const double nevents        = onspill_time/1.695e-6; // N(events) in a week
  const double npot_per_event = 1.6e7*(1.5/3.8); // N(POT) per event
  const double npot           = nevents*npot_per_event; // N(POT) in a week
  const double nmuons         = npot*nmuons_per_pot_run1b_;
  plot_npot_     = npot;
  plot_livetime_ = livetime_week_;
  plot_nmuons_   = nmuons;

  sig_skim_eff_ = datasets["fgam"].nDigi / datasets["fgam"].nGen;
  norm_sig_ = getNorm(datasets["fgam"], f_sig, npot, livetime_week_);
  norm_bkg_ = getNorm(datasets["mnbs"], f_bkg, nevents, livetime_week_);

  // Print summary information
  printf("============================================================\n");
  printf("Livetime      = %.2e s\n", plot_livetime_);
  printf("N(POT)        = %.2e\n"  , plot_npot_);
  printf("N(muon stops) = %.2e\n"  , plot_nmuons_);
  printf("N(events)     = %.2e\n"  , nevents);
  printf("============================================================\n");

  // Set the list of processes to consider
  processes_.clear();
  processes_ = {
    {"RMC"                 , datasets["fgam"], f_sig, norm_sig_,                                0, true , kBlue},
    {"RMC_pu"              , datasets["fgam"], f_sig, norm_sig_,                              100, true , kBlue},
    {"RMC_cpu"             , datasets["fgam"], f_sig, norm_sig_,                              200, true , kBlue},
    {"Cosmics"             , datasets["csms"], f_csm, getNorm(datasets["csms"], f_csm, npot, livetime_week_),   0, false, kGreen-6},
    {"Neutrons"            , datasets["neut"], f_neu, getNorm(datasets["neut"], f_neu, npot, livetime_week_),   0, false, kViolet+6},
    {"Low pileup clusters" , datasets["mnbs"], f_bkg, norm_bkg_,                                0, false, kPink},
    {"Other pileup"        , datasets["mnbs"], f_bkg, norm_bkg_,                              100, false, kViolet},
    {"Calo muon stops"     , datasets["mnbs"], f_bkg, norm_bkg_,                              200, false, kOrange}
  };

  printf("%25s %10s %10s %10s %10s %10s %15s %10s\n", "Process", "N(sampled)", "N(digi)", "N(gen)", "Bare norm", "Norm", "Dataset", "Set offset");
  for(const auto& process : processes_) {
    printf("%25s %10.2e %10.2e %10.2e%10.2e %10.2e %15s %10d\n",
           process.name.Data(),
           getNSampled(process.f),
           process.dataset.nDigi,
           process.dataset.nGen,
           process.dataset.norm((process.dataset.name.BeginsWith("mnbs")) ? nevents : plot_npot_, plot_livetime_),
           process.norm,
           process.dataset.name.Data(),
           process.set_offset);
  }
  // Set up the figure directory and style
  dir_ = (tag) ? Form("figures/rmc_vs_bkg_nt_%s", tag) : "figures/rmc_vs_bkg";
  gSystem->Exec(Form("mkdir -p %s", dir_.Data()));
  gStyle->SetOptStat(0);

  const double emin = 60.;
  const double emax = 120.;

  // Plot by process
  // Plot the histograms
  vector<int> proc_sets = {70, 71, 72, 73, 74};
  for(const int set : proc_sets) {
    for(const bool normalize : {false, true}) {
      plot("cluster_energy"                 , set, normalize, 2, emin,  emax, true, false);
      plot("cluster_time"                   , set, normalize, 5, 600., 1650.);
      plot("cluster_radius"                 , set, normalize, 1, 300.,  700.);
      plot("cluster_disk"                   , set, normalize, 1,   0.,    2.);
      plot("cluster_frac_1"                 , set, normalize, 1,   1.,   -1.);
      plot("cluster_frac_2"                 , set, normalize, 1,   1.,   -1.);
      plot("cluster_e9_over_e"              , set, normalize, 1,   0.,    1.);
      plot("cluster_e25_over_e"             , set, normalize, 1,   0.,    1.);
      plot("cluster_e8_over_e"              , set, normalize, 1,   0.,    1.);
      plot("cluster_e24_over_e"             , set, normalize, 1,   0.,    1.);
      plot("cluster_ncr"                    , set, normalize, 1,   0.,   15.);
      plot("cluster_second_moment"          , set, normalize, 5,   1.,   -1.);
      plot("cluster_t_var"                  , set, normalize, 2,   0.,    5.);
      plot("time_cluster_nhits"             , set, normalize, 1,   0.,  100.);
      plot("time_cluster_nstraw_hits"       , set, normalize, 1,   0.,  100.);
      plot("time_cluster_nhigh_z_hits"      , set, normalize, 1,   0.,   20.);
      plot("crv_dt_corrected"               , set, normalize, 2, -100.,  100.);
      plot("sim_1_2_nhits"                  , set, normalize, 1,   1.,   -1.);
      plot("sim_1_edep"                     , set, normalize, 1,   0.,  100.);
      plot("sim_2_edep"                     , set, normalize, 1,   0.,  100.);
      plot("sim_1_type"                     , set, normalize, 1,  -1.,   10.);
    }
  }

  // Plot the histograms
  vector<int> sets = {0, 1, 3, 4, 18, 20, 23, 24, 35, 38};
  for(const int set : sets) {
    plot_gen_eff(f_sig, set);
    plot_signal(f_sig, "cluster_energy", set, 1,  emin, emax);
    plot_signal(f_sig, "cluster_time"  , set, 5, 600., 1650.);
    for(const bool normalize : {false, true}) {
      plot("cluster_energy"                 , set, normalize, 2, emin,  emax, f_sig, f_bkg);
      plot("cluster_time"                   , set, normalize, 5, 600., 1650., f_sig, f_bkg);
      plot("cluster_radius"                 , set, normalize, 1, 300.,  700., f_sig, f_bkg);
      plot("cluster_disk"                   , set, normalize, 1,   0.,    2., f_sig, f_bkg);
      plot("cluster_frac_1"                 , set, normalize, 1,   1.,   -1., f_sig, f_bkg);
      plot("cluster_frac_2"                 , set, normalize, 1,   1.,   -1., f_sig, f_bkg);
      plot("cluster_ncr"                    , set, normalize, 1,   0.,   15., f_sig, f_bkg);
      plot("cluster_second_moment"          , set, normalize, 5,   1.,   -1., f_sig, f_bkg);
      plot("cluster_t_var"                  , set, normalize, 2,   0.,    5., f_sig, f_bkg);
      plot("time_cluster_nhits"             , set, normalize, 1,   0.,  100., f_sig, f_bkg);
      plot("time_cluster_nstraw_hits"       , set, normalize, 1,   0.,  100., f_sig, f_bkg);
      plot("time_cluster_nhigh_z_hits"      , set, normalize, 1,   0.,   20., f_sig, f_bkg);
      plot("sim_1_2_nhits"                  , set, normalize, 1,   1.,   -1., f_sig, f_bkg);
      plot("sim_1_edep"                     , set, normalize, 1,   0.,  100., f_sig, f_bkg);
      plot("sim_2_edep"                     , set, normalize, 1,   0.,  100., f_sig, f_bkg);
    }
  }

}
