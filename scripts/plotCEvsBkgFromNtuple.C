// Plot CE vs. Bkg

#include "Run1BAna/scripts/plotSigvsBkgFromNtuple.C"

double br_sig_;

//------------------------------------------------------------------------------
void plotCEvsBkgFromNtuple(const char* tag = "v04") {

  // Open the data files
  TFile* f_sig = TFile::Open("Run1BAna.cele4b1s51r0001.hist", "READ");
  TFile* f_bkg = TFile::Open("Run1BAna.mnbs4b1s51r0001.hist", "READ");
  TFile* f_dio = TFile::Open("Run1BAna.diob4b1s51r0000.hist", "READ");
  if (!f_sig || f_sig->IsZombie() ||
      !f_bkg || f_bkg->IsZombie() ||
      !f_dio || f_dio->IsZombie()
      ) {
    Error(__func__, "Could not open files!");
    return;
  }

  // Get the normalization information
  nnt_sig_ = getNSampled(f_sig, -1); // N(events) in the input dataset
  nnt_bkg_ = getNSampled(f_bkg, -1);
  const int nnt_dio = getNSampled(f_dio, -1);
  const int n_saved_sig = getNSampled(f_sig, 0); // N(events) in the ntuples
  const int n_saved_bkg = getNSampled(f_bkg, 0);
  const int n_saved_dio = getNSampled(f_dio, 0);
  if(nnt_sig_ <= 0. || nnt_bkg_ <= 0. || nnt_dio <= 0.) {
    Error(__func__, "Invalid normalization!");
    return;
  }
  br_sig_  = 1.e-9*muon_capture_fraction_; // signal branching fraction
  const double nevents = livetime_week_*duty_cycle_1bb_/1.695e-6; // N(events) in a week
  const double npot    = nevents*1.6e7*(1.5/3.8); // 0.375;
  const double nmuons  = npot*nmuons_per_pot_run1b_;
  const double nsignal = (1.-muon_capture_fraction_)*br_sig_*nmuons;
  sig_skim_eff_ = (1691071./2000000000); // digi dataset
  // sig_skim_eff_ = (); // mcs dataset
  const double dio_skim_eff = (11991587. / 2000000000.); // digi dataset
  const double norm_s  = 1. * sig_skim_eff_; // sample creation + filtering factors
  const double norm_b = 1.; //nnt_bkg_ / getNSampled(f_bkg, -1); // mcs cluster skim
  const double ndio = dio_frac_60_80_ * (1. - muon_capture_fraction_) * nmuons; // N(DIO in 60-80 MeV)
  norm_sig_ = nsignal*norm_s/nnt_sig_;
  norm_bkg_ = nevents*norm_b/nnt_bkg_;
  const double sig_nt_eff = n_saved_sig / nnt_sig_;
  const double ses = 1./(nmuons*muon_capture_fraction_*sig_skim_eff_* sig_nt_eff); // single event sensitivity
  const double norm_dio = ndio*dio_skim_eff/nnt_dio;
  std::cout << "N(sampled): CE = " << nnt_sig_ << " Bkg = " << nnt_bkg_ << " DIO = " << nnt_dio << std::endl;
  std::cout << "N(saved): CE = " << n_saved_sig << " Bkg = " << n_saved_bkg << " DIO = " << n_saved_dio << std::endl;
  std::cout << "Eff(dataset): CE = " << sig_skim_eff_ << " Bkg = " << norm_b << " DIO = " << dio_skim_eff << std::endl;
  std::cout << "N(POT) = " << npot << " N(muons) = " << nmuons << " N(CE) = " << nsignal << " N(DIO) = " << ndio << std::endl;
  std::cout << "Norms: CE = " << norm_sig_ << " Bkg = " << norm_bkg_ << " DIO = " << norm_dio << std::endl;
  std::cout << "N(CE | 0) = " << getNSampled(f_sig, 0)*norm_sig_ << " N(Bkg | 0) = " << getNSampled(f_bkg, 0)*norm_bkg_ << std::endl;
  std::cout << "I(CE | 0) = " << getIntegral(f_sig, 0)*norm_sig_ << " I(Bkg | 0) = " << getIntegral(f_bkg, 0)*norm_bkg_ << std::endl;
  std::cout << "Eff(CE ntuple) = " << sig_nt_eff << std::endl;
  std::cout << "SES = " << ses << std::endl;

  // Set the list of processes to consider
  processes_ = {
    {"#mu^{-}#rightarrowe^{-}", f_sig, norm_sig_,   0, true , kBlue},
    {"CE_pu"   , f_sig, norm_sig_, 100, true , kBlue},
    {"CE_cpu"  , f_sig, norm_sig_, 200, true , kBlue},
    {"DIO tail", f_dio, norm_dio ,   0, false, kGreen-6},
    {"DIO PU"  , f_bkg, norm_bkg_,   0, false, kPink},
    {"Pileup"  , f_bkg, norm_bkg_, 100, false, kViolet},
    {"CaloMu"  , f_bkg, norm_bkg_, 200, false, kOrange}
  };

  // Set up the figure directory and style
  dir_ = (tag) ? Form("figures/ce_vs_bkg_nt_%s", tag) : "figures/ce_vs_bkg_nt";
  gSystem->Exec(Form("mkdir -p %s", dir_.Data()));
  gStyle->SetOptStat(0);

  // Plot by process
  vector<int> proc_sets = {70, 71};
  for(const int set : proc_sets) {
    for(const bool normalize : {false, true}) {
      plot("cluster_energy"                 , set, normalize, 2,  70.,  120.);
      plot("cluster_time"                   , set, normalize, 5, 400., 1800.);
      plot("cluster_radius"                 , set, normalize, 1, 300.,  700.);
      plot("cluster_disk"                   , set, normalize, 1,   0.,    2.);
      plot("cluster_frac_1"                 , set, normalize, 1,   1.,   -1.);
      plot("cluster_frac_2"                 , set, normalize, 1,   1.,   -1.);
      plot("cluster_ncr"                    , set, normalize, 1,   0.,   10.);
      plot("cluster_second_moment"          , set, normalize, 5,   1.,   -1.);
      plot("cluster_t_var"                  , set, normalize, 1,   0.,    5.);
      plot("time_cluster_nhits"             , set, normalize, 2,   0.,  100.);
      plot("time_cluster_nstraw_hits"       , set, normalize, 2,   0.,  100.);
      plot("time_cluster_nhigh_z_hits"      , set, normalize, 1,   0.,   20.);
      plot("sim_1_2_nhits"                  , set, normalize, 1,   1.,   -1.);
      plot("sim_1_edep_frac"                , set, normalize, 1,   0.,    1.);
      plot("sim_2_edep_frac"                , set, normalize, 1,   0.,    1.);
      plot("sim_1_type"                     , set, normalize, 1,  -1.,   10.);
      plot("sim_1_pdg"                      , set, normalize, 1, -15.,   15.);
    }
  }

  // Plot the histograms
  vector<int> sets = {0, 3, 20, 23, 30};
  for(const int set : sets) {
    plot_gen_eff(f_sig, set);
    plot_signal(f_sig, "cluster_energy", set, 1,  60.,  110.);
    plot_signal(f_sig, "cluster_time"  , set, 5, 450., 1650.);
    for(const bool normalize : {false, true}) {
      plot("cluster_energy"                 , set, normalize, 1,  60.,  120., f_sig, f_bkg);
      plot("cluster_time"                   , set, normalize, 2, 400., 1650., f_sig, f_bkg);
      plot("cluster_radius"                 , set, normalize, 1, 300.,  700., f_sig, f_bkg);
      plot("cluster_disk"                   , set, normalize, 1,   0.,    2., f_sig, f_bkg);
      plot("cluster_frac_1"                 , set, normalize, 1,   1.,   -1., f_sig, f_bkg);
      plot("cluster_frac_2"                 , set, normalize, 1,   1.,   -1., f_sig, f_bkg);
      plot("cluster_ncr"                    , set, normalize, 1,   0.,   10., f_sig, f_bkg);
      plot("cluster_second_moment"          , set, normalize, 5,   1.,   -1., f_sig, f_bkg);
      plot("cluster_t_var"                  , set, normalize, 1,   0.,    5., f_sig, f_bkg);
      plot("time_cluster_nhits"             , set, normalize, 2,   0.,  100., f_sig, f_bkg);
      plot("time_cluster_nstraw_hits"       , set, normalize, 2,   0.,  100., f_sig, f_bkg);
      plot("time_cluster_nhigh_z_hits"      , set, normalize, 1,   0.,   20., f_sig, f_bkg);
      plot("sim_1_2_nhits"                  , set, normalize, 1,   1.,   -1., f_sig, f_bkg);
      plot("sim_1_edep_frac"                , set, normalize, 1,   0.,    1., f_sig, f_bkg);
      plot("sim_2_edep_frac"                , set, normalize, 1,   0.,    1., f_sig, f_bkg);
    }
  }

}
