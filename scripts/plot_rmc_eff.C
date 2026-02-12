// Plot efficiencies for the flat gamma sample


//-----------------------------------------------------------------------------------------------
TH1* getHistogram(TFile* f, const char* dir, const char* name) {
  TH1* h = (TH1*) f->Get(Form("%s/%s", dir, name));
  if(h) h = (TH1*) h->Clone(Form("%s_%s", name, dir));
  return h;
}

//-----------------------------------------------------------------------------------------------
void printPlot(TH1* h, const char* name,
               const char* outdir,
               int rebin = 1, double xmin = 1., double xmax = -1.,
               const char* title=nullptr) {

  TCanvas c(name, name, 800, 600);
  if(rebin > 1) h->Rebin(rebin);
  h->Draw("hist");
  h->SetLineColor(kBlue);
  h->SetLineWidth(3);
  // h->SetFillColor(kAtlantic);
  h->SetFillColor(kBlue);
  h->SetFillStyle(3004);
  if(xmin < xmax) h->GetXaxis()->SetRangeUser(xmin, xmax);
  if(title) h->SetTitle(title);

  gPad->Update();
  c.SaveAs(Form("%s/%s.png", outdir, name));
  c.SetLogy();
  c.SaveAs(Form("%s/%s_log.png", outdir, name));
}

//-----------------------------------------------------------------------------------------------
void printPlot(TFile* f,
               const char* dir, const int set, const char* name,
               const char* outdir,
               int rebin = 1, double xmin = 1., double xmax = -1.,
               const char* title=nullptr) {
  TH1* h  = getHistogram(f, Form("Run1BAna/%s_%i", dir, set), name);
  if(!h) {
    std::cerr << "Histogram " << dir << "/" << name << " not found.\n";
    return;
  }
  printPlot(h, Form("%s_%s_%i", name, dir, set), outdir, rebin, xmin, xmax, title);
}

//-----------------------------------------------------------------------------------------------
void print_ratio_plot(TFile* f,
                      const char* dir_1, const char* name_1,
                      const char* dir_2, const char* name_2,
                      const char* outdir, const char* outname,
                      int rebin = 1, double xmin = 1., double xmax = -1.,
                      const char* title=nullptr) {
  TH1* h_1  = getHistogram(f, dir_1, name_1);
  TH1* h_2  = getHistogram(f, dir_2, name_2);
  if(!h_1) {
    std::cerr << "Histogram " << dir_1 << "/" << name_1 << " not found.\n";
    return;
  }
  if(!h_2) {
    std::cerr << "Histogram " << dir_2 << "/" << name_2 << " not found.\n";
    return;
  }
  TH1* h = (TH1*) h_1->Clone(Form("%s_ratio", h_1->GetName()));
  h->Divide(h_2);
  printPlot(h, outname, outdir, rebin, xmin, xmax, title);
}

//-----------------------------------------------------------------------------------------------
void plot_energy_ratio_vs_trig_path(TFile* f,
                                     int set,
                                     const char* outdir,
                                     int rebin = 1) {
  // Retrieve 2D histograms energy_vs_trig_path for both sets
  TH2* h_2d = (TH2*) f->Get(Form("Run1BAna/sim_%i/energy_vs_trig_path", set));

  if(!h_2d) {
    std::cerr << "Cannot find energy_vs_trig_path histogram for set " << set << "\n";
    return;
  }
  std::string dir = Form("%s/trig_%i", outdir, set);
  gSystem->Exec(TString::Format("mkdir -p %s", dir.c_str()));

  int nTrigPaths = h_2d->GetNbinsX();

  // Use the Offline trigger path to normalize the efficiency for each trigger path
  TH1* h_den = nullptr;
  for(int i = 1; i <= nTrigPaths; ++i) {
    TString trigName = h_2d->GetXaxis()->GetBinLabel(i);
    if(trigName == "RecoPath") {
      h_den = h_2d->ProjectionY("h_den", i, i);
      break;
    }
  }
  if(!h_den) {
    std::cerr << "Could not find 'RecoPath' in energy_vs_trig_path histogram for set " << set << "\n";
    return;
  }
  if(h_den->GetEntries() == 0) {
    std::cerr << "'RecoPath' bin in energy_vs_trig_path histogram for set " << set << " is empty.\n";
    return;
  }

  // For each trigger path bin, create a 1D energy ratio plot
  for(int iTrig = 1; iTrig <= nTrigPaths; ++iTrig) {
    TString trigName = h_2d->GetXaxis()->GetBinLabel(iTrig);
    if(trigName == "" || !trigName.Contains("_")) continue;

    // Project Y axis (energy) for this trigger path bin
    TH1* h_num = h_2d->ProjectionY(Form("py_num_%i", iTrig), iTrig, iTrig);

    // Skip if either histogram is empty
    if(h_num->GetEntries() == 0 ) {
      delete h_num;
      continue;
    }

    // Create ratio histogram
    TH1* ratio = (TH1*) h_num->Clone(Form("ratio_trig%i", iTrig));
    ratio->Divide(h_den);
    ratio->GetYaxis()->SetRangeUser(1.e-4, 1.2); // Set Y range for ratio plot
    ratio->SetTitle(Form("Trigger efficiency for path: %s;Generated energy (MeV);Efficiency / %.1g MeV",
                          trigName.Data(), ratio->GetXaxis()->GetBinWidth(1)));

    // Save the ratio plot
    std::string outname = Form("energy_ratio_%s", trigName.Data());
    printPlot(ratio, outname.c_str(), dir.c_str(), rebin, 70., 110.);

    // Clean up
    delete h_num;
    delete ratio;
  }

  std::cout << "Trigger path energy ratio plots for sim_" << set
            << " saved to: " << dir << "\n";
}

//-----------------------------------------------------------------------------------------------
void plot_gen_eff(TFile* f, int set, const double scale, const char* outdir) {
  TH1* h_npot = getHistogram(f, "Run1BAna/evt_0", "npot"); // get N(events) processed
  if(!h_npot) {
    std::cerr << "Cannot find npot histogram for set " << set << "\n";
    return;
  }
  const int nentries = h_npot->GetEntries();
  TH1* h_gen = getHistogram(f, Form("Run1BAna/sim_%i", set), "energy_start"); // get generated energy distribution
  if(!h_gen) {
    std::cerr << "Cannot find energy_start histogram for set " << set << "\n";
    return;
  }
  // Scale the generated energy histogram by the number of events and the provided scale factor
  h_gen->Scale(scale / nentries);
  h_gen->Rebin(2); // Rebin to reduce fluctuations

  // Assume generation was flat between 50 and 110 MeV
  h_gen->Scale((110. - 50.) / h_gen->GetXaxis()->GetBinWidth(1));
  printPlot(h_gen, Form("gen_energy_eff_%i", set), outdir, 1, 70., 110., Form("Generated energy distribution;Energy (MeV);Efficiency / %.1g MeV", h_gen->GetXaxis()->GetBinWidth(1)));
}

//-----------------------------------------------------------------------------------------------
void plot_rmc_eff(const char* infile = "nts.mu2e.FlatGammaMixLowTriggerable.Run1Bab2_best_v1_2.root", const char* outdir = "figures/rmc/eff") {
  gSystem->Exec(TString::Format("mkdir -p %s", outdir));
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(infile, "READ");
  if(!f || f->IsZombie()) {
    std::cerr << "Cannot open file: " << infile << "\n";
    return;
  }

  // scale for nts.mu2e.FlatGammaMixLowTriggerable.Run1Bab2_best_v1_2.root
  const double gen_to_digi_scale = 3189865./969500000.; // N(digi) / N(gen)

  vector<int> event_sets = {0, 10, 20};
  for(int set : event_sets) {
    printPlot(f, "evt", set, "nclusters"   , outdir, 1,     0.,   30.); //, "N(calo clusters)"         );
    printPlot(f, "evt", set, "npot"        , outdir, 1,     1.,   -1.); //, "N(POT)"                   );
    printPlot(f, "sim", set, "energy_start", outdir, 1,    70.,  110.); //, "Generated energy"         );
  }

  vector<int> cluster_sets = {0, 1, 7, 10, 20};
  for(int set : cluster_sets) {
    printPlot(f, "cls", set, "energy", outdir, 1,  50.,  120.); //, "Cluster energy;Energy (MeV)");
    printPlot(f, "cls", set, "time"  , outdir, 1, 300., 1800.); //, "Cluster time;Time (ns)");
    // Absolute
    plot_gen_eff(f, set, gen_to_digi_scale, outdir);
  }

  // Trigger efficiencies
  plot_energy_ratio_vs_trig_path(f, 20, outdir, 1);
}
