// Make standard plots
#include "Run1BAna/analysis/Plotter.C"
Plotter* plotter_ = nullptr;

//---------------------------------------------------------------------------------------------------
// Print the relevant process codes in the model
int print_proc_info(const int selection) {
  if(!plotter_) return 1;
  TString hist = "primary_code";
  TString type = "evt";
  auto backgrounds = plotter_->get_histograms(hist, type, selection,  1);
  auto signals     = plotter_->get_histograms(hist, type, selection, -1);
  auto datas       = plotter_->get_histograms(hist, type, selection,  0);

  printf("Process information for set %4i:\n", selection);
  for(auto h : backgrounds) {
    cout << "Background " << h->GetTitle() << ":\n";
    for(int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
      if(h->GetBinContent(ibin) > 0.) printf("  Process %3i: rate = %10g\n", (int) (h->GetBinCenter(ibin) + 0.5), h->GetBinContent(ibin));
    }
  }

  for(auto h : signals) {
    cout << "Signal " << h->GetTitle() << ":\n";
    for(int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
      if(h->GetBinContent(ibin) > 0.) printf("  Process %3i: rate = %10g\n", (int) (h->GetBinCenter(ibin) + 0.5), h->GetBinContent(ibin));
    }
  }

  for(auto h : datas) {
    cout << "Data " << h->GetTitle() << ":\n";
    for(int ibin = 1; ibin <= h->GetNbinsX(); ++ibin) {
      if(h->GetBinContent(ibin) > 0.) printf("  Process %3i: rate = %10g\n", (int) (h->GetBinCenter(ibin) + 0.5), h->GetBinContent(ibin));
    }
  }

  return 0;
}

//---------------------------------------------------------------------------------------------------
int make_plots(const bool mumem = true,
               vector<int> sets = {0, 10},
               TString dataset = "", TString tag = "") {
  plotter_ = new Plotter();
  plotter_->bkgs_ = {"dio_90_inf", "dio_80_90", "dio_60_80"}; // backgrounds included in the plotting
  // plotter_->bkgs_ = {"dio_90_inf", "dio_80_90"}; // backgrounds included in the plotting

  TString figdir = (mumem) ? "figures/plots/mumem" : "figures/plots/mumep";
  if(dataset != "") figdir += "_" + dataset;
  if(tag != "") figdir += "_" + tag;
  plotter_->figdir_ = figdir;
  if(plotter_->init(mumem, dataset, tag)) return 1;
  if(dataset.BeginsWith("mds1")) plotter_->stack_signal_ = 1; // include the signal in the stacks
  printf("------------------------------------------------------\n");
  printf("Normalization: N(POT) = %.2e, livetime = %.2e\n", npot_, livetime_);
  printf("------------------------------------------------------\n");

  // return 0;

  const bool mds = dataset.Contains("mds");
  const double e_min(70.), e_max(110.);
  TCanvas* c;
  const double base_br(signal_br_);
  plotter_->update_signal_br(signal_br_);
  plotter_->use_offsets_ = false; //don't use control regions for initial counts
  for(int set : sets) {
    if(set < 0) continue;
    for(int logy = 0; logy < 2; ++logy) {
      c = plotter_->print_stack(plot_t("energy"                 , "cls", set, 2, e_min, e_max, 1., -1., logy, false, "Cluster energy", "MeV"));
      c = plotter_->print_stack(plot_t("time"                   , "cls", set, 1,   0., 2000. , 1., -1., logy, false, "Cluster time", "ns"));
      c = plotter_->print_stack(plot_t("radius"                 , "cls", set, 1, 300.,  700. , 1., -1., logy, false, "Cluster radius", "mm"));
      c = plotter_->print_stack(plot_t("ncr"                    , "cls", set, 0,   0.,   10. , 1., -1., logy, false, "N(crystals)", ""));
      c = plotter_->print_stack(plot_t("disk"                   , "cls", set, 0,   1.,   -1. , 1., -1., logy, false, "Cluster disk", ""));
      c = plotter_->print_stack(plot_t("line_dt"                , "cls", set, 1, -50.,   50. , 1., -1., logy, false, "T(cluster) - T(line)", "ns"));
      c = plotter_->print_stack(plot_t("line_dr"                , "cls", set, 1,   0.,  500. , 1., -1., logy, false, "|cluster(x,y) - line(x,y)|", "mm"));
      c = plotter_->print_stack(plot_t("nmatched_lines"         , "cls", set, 1,   0.,    5. , 1., -1., logy, false, "N(matched lines)", ""));
      c = plotter_->print_stack(plot_t("energy_per_crystal"     , "cls", set, 0,   0.,  100. , 1., -1., logy, false, "E / N(crystals)", "MeV"));
      c = plotter_->print_stack(plot_t("frac_first_crystal"     , "cls", set, 0,   0.,    1. , 1., -1., logy, false, "E_{1} / E", ""));
      c = plotter_->print_stack(plot_t("frac_first_two_crystals", "cls", set, 0,   0.,    1. , 1., -1., logy, false, "(E_{1} + E_{2}) / E", ""));
      c = plotter_->print_stack(plot_t("energy_start"           , "sim", set, 1,  50.,  110. , 1., -1., logy, false, "Cluster energy", "MeV"));
      c = plotter_->print_stack(plot_t("nclusters"              , "evt", set, 0,   0.,   10. , 1., -1., logy, false, "N(clusters)", ""));
      c = plotter_->print_stack(plot_t("ngood_clusters"         , "evt", set, 0,   0.,   10. , 1., -1., logy, false, "N(clusters | ID)", ""));
      c = plotter_->print_stack(plot_t("n_lines"                , "evt", set, 0,   0.,   10. , 1., -1., logy, false, "N(lines)", ""));
      c = plotter_->print_stack(plot_t("ngood_lines"            , "evt", set, 0,   0.,   10. , 1., -1., logy, false, "N(lines | ID)", ""));
    }
  }
  return 0;
}
