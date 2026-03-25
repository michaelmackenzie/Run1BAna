// Plot RPC vs. Bkg

#include "Run1BAna/scripts/plotSigvsBkgFromNtuple.C"

//------------------------------------------------------------------------------
void plotRPCvsBkgFromNtuple(const char* tag = "v05") {

  TString csm_file = "Run1BAna.csms4b0s51r0002.hist";

  // Open the data files
  TFile* f_sig = TFile::Open("Run1BAna.rpce4b0s51r0002.hist", "READ");
  TFile* f_bkg = TFile::Open("Run1BAna.mnbs4b1s51r0002.hist", "READ");
  TFile* f_csm = TFile::Open(csm_file, "READ");
  if (!f_bkg || f_bkg->IsZombie() ||
      !f_sig || f_sig->IsZombie() ||
      !f_csm || f_csm->IsZombie()
      ) {
    Error(__func__, "Could not open files!");
    return;
  }

  // Get the normalization information
  nnt_sig_ = getNSampled(f_sig, -1); // N(events) in the input dataset
  nnt_bkg_ = getNSampled(f_bkg, -1);
  const int nnt_csm = getNSampled(f_csm, -1);
  const int n_saved_sig = getNSampled(f_sig, 0); // N(events) in the ntuples
  const int n_saved_bkg = getNSampled(f_bkg, 0);
  const int n_saved_csm = getNSampled(f_csm, 0);
  if(nnt_sig_ <= 0. || nnt_bkg_ <= 0.) {
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

  // RPC info
  const double rpc_skim_eff = (28757. / 1002530508.) * (9487. / 28757.); // dts dataset * (digi / dts)
  const double rpc_stops = (16096977. /  100000000.); // PiTargetStops eff
  const double rpc_beam  = (11978542. / 1000000000.); // PiBeam eff
  const double nrpc = npot * rpc_beam * rpc_stops * rpc_br_; // N(infinite lifetime RPC)
  const double norm_rpc = nrpc*rpc_skim_eff/nnt_sig_;
  sig_skim_eff_ = rpc_skim_eff;
  norm_sig_ = norm_rpc;

  // Cosmic info
  const double livetime_digi = 883457; // (2377000000. / 2585823777.) * 556000.; // N(gen digi) / N(gen sim) * livetime (sim)
  const double ndigi = 115953402; //55369216.; // N(events) in the digi dataset
  const double norm_csm = onspill_time / livetime_digi * (ndigi / nnt_csm);

  // Pileup info
  const double norm_b = 1.;
  norm_bkg_ = nevents*norm_b/nnt_bkg_;

  // Print info
  std::cout << "N(sampled): RPC = " << nnt_sig_ << " Bkg = " << nnt_bkg_ << std::endl;
  std::cout << "N(saved): RPC = " << n_saved_sig << " Bkg = " << n_saved_bkg << std::endl;
  std::cout << "Eff(dataset): RPC = " << sig_skim_eff_ << " Bkg = " << norm_b << std::endl;
  std::cout << "N(POT) = " << npot << " N(muons) = " << nmuons << " N(RPC) = " << nrpc << std::endl;
  std::cout << "Norms: RPC = " << norm_sig_ << " Bkg = " << norm_bkg_ << std::endl;
  std::cout << "N(RPC | 0) = " << getNSampled(f_sig, 0)*norm_sig_ << " N(Bkg | 0) = " << getNSampled(f_bkg, 0)*norm_bkg_ << std::endl;
  std::cout << "I(RPC | 0) = " << getIntegral(f_sig, 0)*norm_sig_ << " I(Bkg | 0) = " << getIntegral(f_bkg, 0)*norm_bkg_ << std::endl;

  std::cout << "\n----------------------------------------------------\n";
  std::cout << "RPC info:\n";
  std::cout << "N(POT) = " << npot << std::endl;
  std::cout << "Eff(PiBeam) = " << rpc_beam << std::endl;
  std::cout << "Eff(PiStops) = " << rpc_stops << std::endl;
  std::cout << "N(RPC infinite lifetime) = " << nrpc << std::endl;
  std::cout << "Eff(Stop time + HitCalo + digi time) = " << rpc_skim_eff << std::endl;
  std::cout << "Eff(cluster) = " << getNSampled(f_sig,0) * 1./nnt_sig_ << std::endl;
  std::cout << "N(RPC cluster) = " << getNSampled(f_sig, 0) * norm_rpc << std::endl;
  std::cout << "N(RPC cluster | time weights) = " << getIntegral(f_sig, 0) * norm_rpc << std::endl;
  std::cout << "----------------------------------------------------\n\n";

  // Set the list of processes to consider
  processes_ = {
    {"RPC"    , f_sig, norm_sig_,   0, true , kBlue},
    {"RPC_pu" , f_sig, norm_sig_, 100, true , kBlue},
    {"RPC_cpu", f_sig, norm_sig_, 200, true , kBlue},
    {"Cosmics", f_csm, norm_csm ,   0, false, kGreen-6},
    {"Low pileup clusters"  , f_bkg, norm_bkg_,   0, false, kPink},
    {"Other pileup"  , f_bkg, norm_bkg_, 100, false, kViolet},
    {"Calo muon stops"  , f_bkg, norm_bkg_, 200, false, kOrange}
  };

  // Set up the figure directory and style
  dir_ = (tag) ? Form("figures/rpc_vs_bkg_nt_%s", tag) : "figures/rpc_vs_bkg";
  gSystem->Exec(Form("mkdir -p %s", dir_.Data()));
  gStyle->SetOptStat(0);

  // Plot by process
  // Plot the histograms
  vector<int> proc_sets = {90, 91, 92, 93, 94, 95, 97};
  for(const int set : proc_sets) {
    for(const bool normalize : {false, true}) {
      plot("cluster_energy"                 , set, normalize, 2,  60.,  140., true, false);
      plot("cluster_time"                   , set, normalize, 2, 450., 1000., true);
      plot("cluster_radius"                 , set, normalize, 1, 300.,  700.);
      plot("cluster_disk"                   , set, normalize, 1,   0.,    2.);
      plot("cluster_frac_1"                 , set, normalize, 1,   1.,   -1.);
      plot("cluster_frac_2"                 , set, normalize, 1,   1.,   -1.);
      plot("cluster_ncr"                    , set, normalize, 1,   0.,   10.);
      plot("cluster_second_moment"          , set, normalize, 5,   1.,   -1.);
      plot("cluster_t_var"                  , set, normalize, 1,   0.,    5.);
      plot("time_cluster_nhits"             , set, normalize, 1,   0.,  100.);
      plot("time_cluster_nstraw_hits"       , set, normalize, 1,   0.,  100.);
      plot("time_cluster_nhigh_z_hits"      , set, normalize, 1,   0.,   20.);
      plot("sim_1_2_nhits"                  , set, normalize, 1,   1.,   -1.);
      plot("sim_1_edep"                     , set, normalize, 1,   0.,  150.);
      plot("sim_2_edep"                     , set, normalize, 1,   0.,  100.);
      plot("sim_1_time"                     , set, normalize, 1, 300., 2000.);
      plot("sim_2_time"                     , set, normalize, 1, 300., 2000.);
      plot("sim_1_type"                     , set, normalize, 1,  -1.,   10.);
    }
  }

  // Plot the histograms
  vector<int> sets = {0};
  for(const int set : sets) {
    plot_gen_eff(f_sig, set);
    plot_signal(f_sig, "cluster_energy", set, 2,  60.,  140.);
    plot_signal(f_sig, "cluster_time"  , set, 5, 400., 1650.);
    for(const bool normalize : {false, true}) {
      plot("cluster_energy"                 , set, normalize, 2,  60.,  140., f_sig, f_bkg);
      plot("cluster_time"                   , set, normalize, 2, 400., 2000., f_sig, f_bkg);
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
      plot("sim_1_type"                     , set, normalize, 1,  -1.,   10., f_sig, f_bkg);
    }
  }

}
