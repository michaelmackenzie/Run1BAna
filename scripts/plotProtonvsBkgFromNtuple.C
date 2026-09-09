// Plot High KE Protons vs. Bkg

#include "Run1BAna/scripts/plotSigvsBkgFromNtuple.C"
#include "Run1BAna/scripts/build_model.C"


//------------------------------------------------------------------------------
void plotProtonvsBkgFromNtuple(const char* tag = "v40", TString hist_tag = "",
                               const double pot_scale = 1.) {

  auto datasets = getDatasets(tag);
  if(datasets.empty()) {
    Error(__func__, "No datasets found for tag %s", tag);
    return;
  }

  const vector<TString> included_dataset_keys = {"prot", "mnbs", "csms", "neut", "fgam", "pgam"};
  map<TString, TFile*> files;
  if(!openIncludedDatasetFiles(datasets, included_dataset_keys, files, hist_tag, __func__)) return;

  TFile* f_sig = getDatasetFile(files, "prot");
  TFile* f_bkg = getDatasetFile(files, "mnbs");
  if(!f_sig || !f_bkg) {
    Error(__func__, "Missing required proton or background datasets in included keys");
    return;
  }

  // General info
  const double onspill_time   = livetime_week_*duty_cycle_1bb_;
  double       nevents        = onspill_time/1.695e-6; // N(events) in a week
  const double npot_per_event = pot_scale*getNPOT(f_sig); // N(POT) per event, from simulated mean value
  const double npot           = nevents*npot_per_event; // N(POT) in a week
  const double nmuons         = npot*nmuons_per_pot_run1b_;
  const float energy =  8.*1.602176634e-10; // 8 GeV proton KE in joules
  const float power = energy * (npot / livetime_week_)/1000.; // in kW
  plot_npot_     = npot;
  plot_livetime_ = livetime_week_;
  plot_nmuons_   = nmuons;

  // Adjust N(events) given POT scaling
  nevents *= pot_scale;

  const vector<TString> enabled_process_ids = {
    "protons", "protons_pu", "protons_cpu", "cosmics", "neutrons", "poly", "rmc", "pileup_lo", "pileup_ot", "calomu"
  };
  const auto process_specs = selectNominalProcessSpecs(enabled_process_ids, false);

  // Set the list of processes to consider
  processes_ = buildProcesses(datasets, files, included_dataset_keys, process_specs, npot, nevents, livetime_week_, rmue_);

  norm_sig_ = getNorm(datasets["prot"], f_sig, npot, livetime_week_);
  norm_bkg_ = getNorm(datasets["mnbs"], f_bkg, nevents, livetime_week_);

  // Print summary information
  printf("============================================================\n");
  printf("Livetime       = %.3e s\n", plot_livetime_);
  printf("N(POT)         = %.3e\n"  , plot_npot_);
  printf("N(muon stops)  = %.3e\n"  , plot_nmuons_);
  printf("N(events)      = %.3e\n"  , nevents);
  printf("N(POT / event) = %.3e\n"  , npot_per_event);
  printf("Power          = %.3f\n"  , power);
  printf("============================================================\n");

  // Print information about each process
  printf("%25s %10s %10s %10s %10s %10s %15s %10s\n", "Process", "N(sampled)", "N(digi)", "N(gen)", "Bare norm", "Norm", "Dataset", "Set offset");
  for(auto& process : processes_) {
    if(pot_scale != 1. && (process.name == "Other pileup" || process.name == "Calo muon stops")) process.norm *= pot_scale; // goes as rate^2
    if(process.name.Contains("Proton")) process.is_signal = true;
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
  dir_ = (tag) ? Form("figures/proton_vs_bkg_nt_%s", tag) : "figures/proton_vs_bkg_nt";
  if(hist_tag != "") dir_ += "_" + hist_tag;
  if(pot_scale != 1.) dir_ += Form("_pot_scale_%id%02i", int(pot_scale), int(pot_scale*100.) - int(pot_scale)*100);
  gSystem->Exec(Form("mkdir -p %s", dir_.Data()));
  gStyle->SetOptStat(0);

  const double emin = 60.;
  const double emax = 100.;
  stack_sig_ = true;
  signal_color_ = kAtlantic;

  // Plot by process
  vector<int> proc_sets = {70,40,43,44};
  for(const int set : proc_sets) {
    for(const bool normalize : {false}) {
      plot("cluster_energy"                 , set, normalize, 2, emin,  emax, "MeV", true);
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
      plot("line_nhits"                     , set, normalize, 1,   0.,  100.);
      plot("line_cos"                       , set, normalize, 1,  0.8,   1.1);
      plot("sim_1_2_nhits"                  , set, normalize, 1,   1.,   -1.);
      plot("sim_1_edep_frac"                , set, normalize, 1,   0.,    1.);
      plot("sim_2_edep_frac"                , set, normalize, 1,   0.,    1.);
      plot("sim_1_type"                     , set, normalize, 1,  -1.,   10.);
      plot("sim_1_pdg"                      , set, normalize, 1, -15.,   15.);
      plot("gen_energy"                     , set, normalize, 2, emin,  emax, "MeV", true);
      plot("gen_ediff"                      , set, normalize, 2, -50.,   10., "MeV", true);
    }
    evaluate_summary(set);
    plotModel(processes_, "cluster_energy", set);
  }

}
