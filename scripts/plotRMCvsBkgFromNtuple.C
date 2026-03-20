// Plot RMC vs. Bkg

#include "Run1BAna/scripts/plotSigvsBkgFromNtuple.C"

//------------------------------------------------------------------------------
void plotRMCvsBkgFromNtuple(const char* filename_sig = "Run1BAna.fgam4b1s51r0000.hist",
                            const char* filename_bkg = "Run1BAna.mnbs4b1s51r0000.hist",
                            const char* tag = "v04") {

  // Open the data files
  TFile* f_sig = TFile::Open(filename_sig, "READ");
  TFile* f_bkg = TFile::Open(filename_bkg, "READ");
  TFile* f_dio = TFile::Open("Run1BAna.diob4b1s51r0000.hist", "READ");
  TFile* f_rpc = TFile::Open("Run1BAna.rpce4b1s51r0001.hist", "READ");
  TFile* f_csm = TFile::Open("Run1BAna.csms4b0s51r0001.hist", "READ");
  if (!f_sig || f_sig->IsZombie() ||
      !f_bkg || f_bkg->IsZombie() ||
      !f_dio || f_dio->IsZombie() ||
      !f_rpc || f_rpc->IsZombie() ||
      !f_csm || f_csm->IsZombie()
      ) {
    Error(__func__, "Could not open files!");
    return;
  }

  // Get the normalization information
  nnt_sig_ = getNSampled(f_sig, -1); // N(events) in the input dataset
  nnt_bkg_ = getNSampled(f_bkg, -1);
  const int nnt_dio = getNSampled(f_dio, -1);
  const int nnt_rpc = getNSampled(f_rpc, -1);
  const int nnt_csm = getNSampled(f_csm, -1);
  const int n_saved_sig = getNSampled(f_sig, 0); // N(events) in the ntuples
  const int n_saved_bkg = getNSampled(f_bkg, 0);
  const int n_saved_dio = getNSampled(f_dio, 0);
  const int n_saved_rpc = getNSampled(f_rpc, 0);
  const int n_saved_csm = getNSampled(f_csm, 0);
  if(nnt_sig_ <= 0. || nnt_bkg_ <= 0. || nnt_dio <= 0. || nnt_rpc <= 0. || nnt_csm <= 0.) {
    Error(__func__, "Invalid normalization!");
    return;
  }

  // General info
  const double onspill_time = livetime_week_*duty_cycle_1bb_;
  const double nevents      = onspill_time/1.695e-6; // N(events) in a week
  const double npot         = nevents*1.6e7*(1.5/3.8);
  const double nmuons       = npot*nmuons_per_pot_run1b_;

  // RMC info
  const double nrmc    = nmuons*muon_capture_fraction_*br_rmc_/rmc_frac_57_; // N(RMC) assuming closure
  sig_skim_eff_ = (1263859. / 1949000000.); // digi dataset
  const double norm_rmc  = (110. - 50.)/90.1 * sig_skim_eff_; // sample creation + filtering factors
  norm_sig_ = nrmc*norm_rmc/nnt_sig_;

  // DIO info
  const double dio_skim_eff = (11991587. / 2000000000.); // digi dataset
  const double ndio = dio_frac_60_80_ * (1. - muon_capture_fraction_) * nmuons; // N(DIO in 60-80 MeV)
  const double norm_dio = ndio*dio_skim_eff/nnt_dio;

  // RPC info
  const double rpc_skim_eff = (1645309. / 519723527.); // digi dataset
  const double rpc_stops = (16096977. /  100000000.); // PiTargetStops eff
  const double rpc_beam  = (11978542. / 1000000000.); // PiBeam eff
  const double nrpc = npot * rpc_beam * rpc_stops * rpc_br_; // N(infinite lifetime RPC)
  const double norm_rpc = nrpc*rpc_skim_eff/nnt_rpc;

  // Cosmic info
  const double livetime_digi = (2377000000. / 2585823777.) * 556000.; // N(gen digi) / N(gen sim) * livetime (sim)
  const double ndigi = 55369216.; // N(events) in the digi dataset
  const double norm_csm = onspill_time / livetime_digi * (ndigi / nnt_csm);

  // Pileup info
  // const double norm_b  = 3880./100000.; // digi trigger cluster skim
  const double norm_b = 1.; //nnt_bkg_ / getNSampled(f_bkg, -1); // mcs cluster skim
  norm_bkg_ = nevents*norm_b/nnt_bkg_;

  // Print the normalization info
  std::cout << "N(POT) = " << npot << " N(muons) = " << nmuons << " On-Spill time = " << onspill_time << std::endl;
  std::cout << "N(RMC) = " << nrmc << " N(DIO) = " << ndio << " N(RPC) = " << nrpc << std::endl;
  std::cout << "N(sampled): RMC = " << nnt_sig_ << " Bkg = " << nnt_bkg_ << " DIO = " << nnt_dio << " RPC = " << nnt_rpc << std::endl;
  std::cout << "N(saved): RMC = " << n_saved_sig << " Bkg = " << n_saved_bkg << " DIO = " << n_saved_dio << " RPC = " << n_saved_rpc << std::endl;
  std::cout << "Eff(dataset): RMC = " << sig_skim_eff_ << " Bkg = " << norm_b << " DIO = " << dio_skim_eff << " RPC = " << rpc_skim_eff << std::endl;
  std::cout << "Norms: RMC = " << norm_sig_ << " Bkg = " << norm_bkg_ << " DIO = " << norm_dio << " RPC = " << norm_rpc << std::endl;
  std::cout << "N(RMC | 0) = " << getNSampled(f_sig, 0)*norm_sig_ << " N(Bkg | 0) = " << getNSampled(f_bkg, 0)*norm_bkg_ << " N(RPC | 0) = " << getNSampled(f_rpc, 0)*norm_rpc << std::endl;
  std::cout << "I(RMC | 0) = " << getIntegral(f_sig, 0)*norm_sig_ << " I(Bkg | 0) = " << getIntegral(f_bkg, 0)*norm_bkg_ << " I(RPC | 0) = " << getIntegral(f_rpc, 0)*norm_rpc << std::endl;

  // Set the list of processes to consider
  processes_ = {
    {"RMC"    , f_sig, norm_sig_,   0, true , kBlue},
    {"RMC_pu" , f_sig, norm_sig_, 100, true , kBlue},
    {"RMC_cpu", f_sig, norm_sig_, 200, true , kBlue},
    {"Cosmics", f_csm, norm_csm ,   0, false, kGreen-6},
    // {"DIO"    , f_dio, norm_dio ,   0, false, kGreen-6},
    {"DIO-pu" , f_bkg, norm_bkg_,   0, false, kPink},
    {"Pileup" , f_bkg, norm_bkg_, 100, false, kViolet},
    {"CaloMu" , f_bkg, norm_bkg_, 200, false, kOrange}
  };

  // Set up the figure directory and style
  dir_ = (tag) ? Form("figures/rmc_vs_bkg_nt_%s", tag) : "figures/rmc_vs_bkg";
  gSystem->Exec(Form("mkdir -p %s", dir_.Data()));
  gStyle->SetOptStat(0);

  // Plot by process
  // Plot the histograms
  vector<int> proc_sets = {70, 71, 72, 73, 74};
  for(const int set : proc_sets) {
    for(const bool normalize : {false, true}) {
      plot("cluster_energy"                 , set, normalize, 1,  70.,  140.);
      plot("cluster_time"                   , set, normalize, 5, 600., 1650.);
      plot("cluster_radius"                 , set, normalize, 1, 300.,  700.);
      plot("cluster_disk"                   , set, normalize, 1,   0.,    2.);
      plot("cluster_frac_1"                 , set, normalize, 1,   1.,   -1.);
      plot("cluster_frac_2"                 , set, normalize, 1,   1.,   -1.);
      plot("cluster_e9_over_e"              , set, normalize, 1,   0.,    1.);
      plot("cluster_e25_over_e"             , set, normalize, 1,   0.,    1.);
      plot("cluster_e8_over_e"              , set, normalize, 1,   0.,    1.);
      plot("cluster_e24_over_e"             , set, normalize, 1,   0.,    1.);
      plot("cluster_ncr"                    , set, normalize, 1,   0.,   10.);
      plot("cluster_second_moment"          , set, normalize, 5,   1.,   -1.);
      plot("cluster_t_var"                  , set, normalize, 1,   0.,    5.);
      plot("time_cluster_nhits"             , set, normalize, 1,   0.,  100.);
      plot("time_cluster_nstraw_hits"       , set, normalize, 1,   0.,  100.);
      plot("time_cluster_nhigh_z_hits"      , set, normalize, 1,   0.,   20.);
      plot("sim_1_2_nhits"                  , set, normalize, 1,   1.,   -1.);
      plot("sim_1_edep"                     , set, normalize, 1,   0.,  100.);
      plot("sim_2_edep"                     , set, normalize, 1,   0.,  100.);
    }
  }

  // Plot the histograms
  vector<int> sets = {0, 1, 3, 4, 18, 20, 23, 24, 35, 38};
  for(const int set : sets) {
    plot_gen_eff(f_sig, set);
    plot_signal(f_sig, "cluster_energy", set, 1,  70.,  110.);
    plot_signal(f_sig, "cluster_time"  , set, 5, 600., 1650.);
    for(const bool normalize : {false, true}) {
      plot("cluster_energy"                 , set, normalize, 1,  70.,  120., f_sig, f_bkg);
      plot("cluster_time"                   , set, normalize, 5, 600., 1650., f_sig, f_bkg);
      plot("cluster_radius"                 , set, normalize, 1, 300.,  700., f_sig, f_bkg);
      plot("cluster_disk"                   , set, normalize, 1,   0.,    2., f_sig, f_bkg);
      plot("cluster_frac_1"                 , set, normalize, 1,   1.,   -1., f_sig, f_bkg);
      plot("cluster_frac_2"                 , set, normalize, 1,   1.,   -1., f_sig, f_bkg);
      plot("cluster_ncr"                    , set, normalize, 1,   0.,   10., f_sig, f_bkg);
      plot("cluster_second_moment"          , set, normalize, 5,   1.,   -1., f_sig, f_bkg);
      plot("cluster_t_var"                  , set, normalize, 1,   0.,    5., f_sig, f_bkg);
      plot("time_cluster_nhits"             , set, normalize, 1,   0.,  100., f_sig, f_bkg);
      plot("time_cluster_nstraw_hits"       , set, normalize, 1,   0.,  100., f_sig, f_bkg);
      plot("time_cluster_nhigh_z_hits"      , set, normalize, 1,   0.,   20., f_sig, f_bkg);
      plot("sim_1_2_nhits"                  , set, normalize, 1,   1.,   -1., f_sig, f_bkg);
      plot("sim_1_edep"                     , set, normalize, 1,   0.,  100., f_sig, f_bkg);
      plot("sim_2_edep"                     , set, normalize, 1,   0.,  100., f_sig, f_bkg);
    }
  }

}
