// Plot RMC vs. Bkg

#include "Run1BAna/scripts/plotSigvsBkgFromNtuple.C"

//------------------------------------------------------------------------------
void plotRMCvsBkgFromNtuple(const char* tag = "v04") {

  TString sig_file = "Run1BAna.fgam4b1s51r0002.hist";
  TString bkg_file = "Run1BAna.mnbs4b1s51r0002.hist";
  TString dio_file = "Run1BAna.diob4b1s51r0002.hist";
  TString rpc_file = "Run1BAna.rpce4b0s51r0002.hist";
  TString csm_file = "Run1BAna.csms4b0s51r0002.hist";
  const bool is_v06 = TString(tag) == "v06"; // 2 cm target + 10 cm poly
  if(is_v06) {
    sig_file = "Run1BAna.fgam6b0s51r0002.hist";
    bkg_file = "Run1BAna.mnbs6b1s51r0002.hist";
    csm_file = "Run1BAna.csms6b0s51r0002.hist";
  }
  TFile* f_sig = TFile::Open(sig_file, "READ");
  TFile* f_bkg = TFile::Open(bkg_file, "READ");
  TFile* f_dio = TFile::Open(dio_file, "READ");
  TFile* f_rpc = TFile::Open(rpc_file, "READ");
  TFile* f_csm = TFile::Open(csm_file, "READ");
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
  plot_npot_     = npot;
  plot_livetime_ = livetime_week_;
  plot_nmuons_   = nmuons;

  // RMC info
  const double nrmc    = nmuons*muon_capture_fraction_*br_rmc_/rmc_frac_57_; // N(RMC) assuming closure
  if(is_v06) {
    sig_skim_eff_ = (1497757. / 2000000000.) * (14892. / 30000.); // (dts / gen) * (digi / dts) dataset
  } else {
    sig_skim_eff_ = (1263859. / 1949000000.); // digi dataset
  }
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
  const double livetime_digi = 883457; // (2377000000. / 2585823777.) * 556000.; // N(gen digi) / N(gen sim) * livetime (sim)
  const double ndigi = 115953402; //55369216.; // N(events) in the digi dataset
  const double norm_csm = onspill_time / livetime_digi * (ndigi / nnt_csm);

  // Pileup info
  const double norm_b = 1.; // no filtering upstream of reco stage
  norm_bkg_ = nevents*norm_b/nnt_bkg_;

  // Print the normalization info
  std::cout << "N(POT) = " << npot << " N(muons) = " << nmuons << " On-Spill time = " << onspill_time << std::endl;
  std::cout << "N(RMC) = " << nrmc << " N(DIO) = " << ndio << " N(RPC) = " << nrpc << std::endl;

  // Print RMC
  std::cout << "-----------------------------------------------\n";
  std::cout << "N(RMC) = " << nrmc << std::endl;
  std::cout << "N(RMC sampled) = " << nnt_sig_ << std::endl;
  std::cout << "N(RMC saved) = " << n_saved_sig << std::endl;
  std::cout << "Eff(RMC dataset) = " << sig_skim_eff_ << std::endl;
  std::cout << "Norm(RMC) = " << norm_sig_ << std::endl;
  std::cout << "N(RMC | clusters) = " << n_saved_sig * norm_sig_ << std::endl;
  std::cout << "I(RMC | clusters) = " << getIntegral(f_sig, 0)*norm_sig_ << std::endl;
  std::cout << "-----------------------------------------------\n";

  // Print Cosmics
  std::cout << "-----------------------------------------------\n";
  std::cout << "Livetime digi = " << livetime_digi << std::endl;
  std::cout << "N(digi) = " << ndigi << std::endl;
  std::cout << "N(Cosmic sampled) = " << nnt_csm << std::endl;
  std::cout << "N(Cosmic saved) = " << n_saved_csm << std::endl;
  std::cout << "Livetime sampled = " << livetime_digi * nnt_csm / ndigi << std::endl;
  std::cout << "Norm(Cosmic) = " << norm_csm << std::endl;
  std::cout << "N(Cosmic | clusters) = " << n_saved_csm * norm_csm << std::endl;
  std::cout << "I(Cosmic | clusters) = " << getIntegral(f_csm, 0)*norm_csm << std::endl;
  std::cout << "-----------------------------------------------\n";

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
    // {"DIO tail", f_dio, norm_dio ,   0, false, kGreen-6},
    {"Low pileup clusters"  , f_bkg, norm_bkg_,   0, false, kPink},
    {"Other pileup"  , f_bkg, norm_bkg_, 100, false, kViolet},
    {"Calo muon stops"  , f_bkg, norm_bkg_, 200, false, kOrange}
  };

  // Set up the figure directory and style
  dir_ = (tag) ? Form("figures/rmc_vs_bkg_nt_%s", tag) : "figures/rmc_vs_bkg";
  gSystem->Exec(Form("mkdir -p %s", dir_.Data()));
  gStyle->SetOptStat(0);

  const double emin = (is_v06) ?  60. :  60.;
  const double emax = (is_v06) ? 120. : 120.;

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
