// Plot CE vs. Bkg

#include "Run1BAna/scripts/plotSigvsBkgFromNtuple.C"


//------------------------------------------------------------------------------
void plotCEvsBkgFromNtuple(const char* tag = "v04") {

  // Open the data files
  TString sig_file = "Run1BAna.cele4b1s51r0001.hist";
  TString bkg_file = "Run1BAna.mnbs4b1s51r0001.hist";
  TString dio_file = "Run1BAna.diob4b1s51r0000.hist";
  TString csm_file = "Run1BAna.csms4b0s51r0001.hist";
  const bool is_v06 = TString(tag) == "v06"; // 2 cm target + 10 cm poly
  if(is_v06) {
    sig_file = "Run1BAna.cele6b1s51r0002.hist";
    bkg_file = "Run1BAna.mnbs6b1s51r0002.hist";
  }
  TFile* f_sig = TFile::Open(sig_file, "READ");
  TFile* f_bkg = TFile::Open(bkg_file, "READ");
  TFile* f_dio = TFile::Open(dio_file, "READ");
  TFile* f_csm = TFile::Open(csm_file, "READ");
  if (!f_sig || f_sig->IsZombie() ||
      !f_bkg || f_bkg->IsZombie() ||
      !f_dio || f_dio->IsZombie() ||
      !f_csm || f_csm->IsZombie()
      ) {
    Error(__func__, "Could not open files!");
    return;
  }

  // Get the normalization information
  nnt_sig_ = getNSampled(f_sig, -1); // N(events) in the input dataset
  nnt_bkg_ = getNSampled(f_bkg, -1);
  const int nnt_dio = getNSampled(f_dio, -1);
  const int nnt_csm = getNSampled(f_csm, -1);
  const int n_saved_sig = getNSampled(f_sig, 0); // N(events) in the ntuples
  const int n_saved_bkg = getNSampled(f_bkg, 0);
  const int n_saved_dio = getNSampled(f_dio, 0);
  const int n_saved_csm = getNSampled(f_csm, 0);
  if(nnt_sig_ <= 0. || nnt_bkg_ <= 0. || nnt_dio <= 0. || nnt_csm <= 0.) {
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

  // CE info
  br_sig_  = 1.e-9*muon_capture_fraction_; // signal branching fraction
  const double nsignal = (1.-muon_capture_fraction_)*br_sig_*nmuons;
  if(is_v06) {
    sig_skim_eff_ = (325238. / 1257000000.); // digi dataset
  } else {
    sig_skim_eff_ = (1691071./2000000000); // digi dataset
  }
  const double norm_s  = 1. * sig_skim_eff_; // sample creation + filtering factors
  norm_sig_ = nsignal*norm_s/nnt_sig_;
  const double sig_nt_eff = n_saved_sig / nnt_sig_;
  const double ses = 1./(nmuons*muon_capture_fraction_*sig_skim_eff_* sig_nt_eff); // single event sensitivity

  // DIO info
  const double dio_skim_eff = (11991587. / 2000000000.); // digi dataset
  const double ndio = dio_frac_60_80_ * (1. - muon_capture_fraction_) * nmuons; // N(DIO in 60-80 MeV)
  const double norm_dio = ndio*dio_skim_eff/nnt_dio;

  // Cosmic info
  // const double livetime_sim  = 556000.; // livetime reported in the sim dataset
  // const double ngen_sim      = 2585823777.; // N(gen events) in the sim dataset
  // const double ngen_digi     = 2377000000.; // N(gen events) in the digi dataset
  // const double livetime_digi = (ngen_digi / ngen_sim) * livetime_sim; // scale the original livetime by the sample factor
  // const double ndigi         = 55369216.; // N(events) in the digi dataset
  const double livetime_digi = 883457; // (2377000000. / 2585823777.) * 556000.; // N(gen digi) / N(gen sim) * livetime (sim)
  const double ndigi = 115953402; //55369216.; // N(events) in the digi dataset
  const double norm_csm      = onspill_time / livetime_digi * (ndigi / nnt_csm);

  // Pileup info
  const double norm_b = 1.; //nnt_bkg_ / getNSampled(f_bkg, -1); // mcs cluster skim
  norm_bkg_ = nevents*norm_b/nnt_bkg_;

  // Print out info
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
  if(is_v06) {
    processes_ = {
      {"#mu^{-}#rightarrowe^{-}", f_sig, norm_sig_,   0, true , kBlue},
      {"CE_pu"   , f_sig, norm_sig_, 100, true , kBlue},
      {"CE_cpu"  , f_sig, norm_sig_, 200, true , kBlue},
      // {"Cosmics" , f_csm, norm_csm ,   0, false, kViolet+6},
      // {"DIO tail", f_dio, norm_dio ,   0, false, kGreen-6},
      {"Low pileup clusters"  , f_bkg, norm_bkg_,   0, false, kPink},
      {"Other pileup"  , f_bkg, norm_bkg_, 100, false, kViolet}
    };
  } else {
    processes_ = {
      {"#mu^{-}#rightarrowe^{-}", f_sig, norm_sig_,   0, true , kBlue},
      {"CE_pu"   , f_sig, norm_sig_, 100, true , kBlue},
      {"CE_cpu"  , f_sig, norm_sig_, 200, true , kBlue},
      // {"Cosmics" , f_csm, norm_csm ,   0, false, kViolet+6},
      // {"DIO tail", f_dio, norm_dio ,   0, false, kGreen-6},
      {"Low pileup clusters"  , f_bkg, norm_bkg_,   0, false, kPink},
      {"Other pileup"  , f_bkg, norm_bkg_, 100, false, kViolet},
      {"Calo muon stops"  , f_bkg, norm_bkg_, 200, false, kOrange}
    };
  }

  // Set up the figure directory and style
  dir_ = (tag) ? Form("figures/ce_vs_bkg_nt_%s", tag) : "figures/ce_vs_bkg_nt";
  gSystem->Exec(Form("mkdir -p %s", dir_.Data()));
  gStyle->SetOptStat(0);

  const double emin = (is_v06) ?  60. :  70.;
  const double emax = (is_v06) ? 100. : 120.;

  // Plot by process
  vector<int> proc_sets = {70, 71, 80};
  for(const int set : proc_sets) {
    for(const bool normalize : {false, true}) {
      plot("cluster_energy"                 , set, normalize, 2, emin,  emax, true);
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
  vector<int> sets = {0, 3, 5, 6, 7, 20, 23, 30};
  for(const int set : sets) {
    plot_signal(f_sig, "cluster_energy", set, 1, emin, emax);
    plot_signal(f_sig, "cluster_time"  , set, 5, 450., 1650.);
    for(const bool normalize : {false, true}) {
      plot("cluster_energy"                 , set, normalize, 2, emin,  emax, f_sig, f_bkg);
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
      plot("line_nhits"                     , set, normalize, 1,   0.,  100., f_sig, f_bkg);
      plot("sim_1_2_nhits"                  , set, normalize, 1,   1.,   -1., f_sig, f_bkg);
      plot("sim_1_edep_frac"                , set, normalize, 1,   0.,    1., f_sig, f_bkg);
      plot("sim_2_edep_frac"                , set, normalize, 1,   0.,    1., f_sig, f_bkg);
    }
  }

}
