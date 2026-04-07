// Evaluate the double edep distribution

//----------------------------------------------------------------------------------------------------------
// Double a histogram by sampling from the original distribution
TH1* mix_histograms(const TH1* h_1, TH1* h_2, TRandom3& rnd) {
  if(!h_1 || !h_2) {
    std::cerr << "Error: Null histogram pointer passed to mix_histograms" << std::endl;
    return nullptr;
  }

  // Create a new histogram with the same binning as the input histograms
  TH1* h_mixed = dynamic_cast<TH1*>(h_1->Clone("mixed_histogram"));
  h_mixed->Reset();

  // Sample each histogram and fill the mixed histogram
  const double bin_width = h_1->GetBinWidth(1); // Assuming uniform binning
  const long max_entries = 1e7;
  for(long i = 0; i < max_entries; ++i) {
    const double x_1 = h_1->GetRandom(&rnd);
    const double x_2 = h_2->GetRandom(&rnd);
    const double doubled = (x_1 + x_2);
    h_mixed->Fill(doubled);
  }
  h_mixed->Scale(1. / max_entries); // Normalize the histogram

  return h_mixed;
}

//----------------------------------------------------------------------------------------------------------
// Main function
int double_edep(vector<TString> file_names,
                vector<double> efficiencies,
                vector<TString> titles,
                const char* run_dir = ".") {

  // Validate the inputs
  const size_t n_bkgs = file_names.size();
  if(n_bkgs == 0) {
    cerr << __func__ << ": No files given!\n";
    return 1;
  }
  if(efficiencies.size() != n_bkgs || titles.size() != n_bkgs) {
    cerr << __func__ << ": File names, efficiencies per POT, and titles must match!\n";
    return 1;
  }

  // Open the histogram files
  vector<TFile*> files;
  for(auto file_name : file_names) {
    TFile* file = TFile::Open(file_name, "READ");
    if(!file || file->IsZombie()) {
      std::cerr << __func__ <<  ": Error opening file: " << file_name << std::endl;
      return 1;
    }
    files.push_back(file);
  }

  // Get the total calo energy deposit per simulated particle histograms
  TH1* h_edep = nullptr;
  vector<TH1*> h_edeps;
  for(size_t index = 0; index < n_bkgs; ++index) {
    TH1* h = dynamic_cast<TH1*>(files[index]->Get("EDepAna/hist_2/total_calo_energy")); // ignore E < 1 MeV steps
    TH1* h_norm = dynamic_cast<TH1*>(files[index]->Get("EDepAna/hist_0/total_calo_energy")); // for absolute rates
    if(!h || !h_norm) {
      std::cerr << "Error: Histogram not found in file: " << files[index]->GetName() << std::endl;
      return 1;
    }
    h = (TH1*) h->Clone(Form("h_edep_%zu", index));
    h->SetDirectory(0);
    h->Scale(efficiencies[index] / h_norm->GetEntries());
    h->SetLineWidth(2);
    if(h_edep) {
      h_edep->Add(h);
    } else {
      h_edep = (TH1*) h->Clone("h_edep");
      h_edep->SetDirectory(0);
    }
    h_edeps.push_back(h);
  }
  for(auto file : files) file->Close();

  // Plot the inputs
  const char* fig_dir = Form("%s/figures");
  gSystem->Exec(Form("mkdir -p %s", fig_dir));
  gStyle->SetOptStat(0);
  TCanvas c("c", "Double Edep", 800, 600);
  TLegend leg_edeps(0.1, 0.8, 0.9, 0.9);
  leg_edeps.SetNColumns(n_bkgs);
  for(size_t index = 0; index < n_bkgs; ++index) {
    if(index == 0) h_edeps[index]->Draw("hist");
    else           h_edeps[index]->Draw("hist same");
    h_edeps[index]->SetLineColor(index + 2);
    leg_edeps.AddEntry(h_edeps[index], titles[index].Data(), "l");
  }
  leg_edeps.Draw();
  c.SaveAs(Form("%s/edeps.png", fig_dir));
  c.SetLogy();
  c.SaveAs(Form("%s/edeps_log.png", fig_dir));
  c.SetLogy(false);

  const double eff   = h_edep->Integral(); // total efficiency per POT
  const double npot  = 1.6e7*(1.5/3.8); // 1.5 kW intensity
  const double n_exp = eff * npot; // expected N(deposits) / average intensity event
  h_edep->Scale(1./eff); // Normalize to probability distribution

  // FIXME: Guess of p(time overlap) and p(space overlap)
  const double p_time_overlap = 0.05;
  const double p_space_overlap = 50.*50. / (3.14*(700.*700. - 350.*350.)); // area of 50x50 cm^2 over area of annulus between 350 and 700 cm radius
  const double p_overlap = p_time_overlap * p_space_overlap;


  // Efficiencies for double and triple edep assuming independent events
  const double n_exp_double = n_exp * p_overlap;
  const double n_exp_triple = n_exp * p_overlap * n_exp_double;
  std::cout << "N(POT) / event average: " << npot << std::endl;
  std::cout << "Efficiency for single edep / POT: " << eff << std::endl;
  std::cout << "Expected per event for single edep: " << n_exp << std::endl;
  std::cout << "Expected per event for double edep: " << n_exp_double << std::endl;
  std::cout << "Expected per event for triple edep: " << n_exp_triple << std::endl;

  // Get the double and triple edep histograms
  TRandom3 rnd(90);
  TH1* h_double = mix_histograms(h_edep  , h_edep, rnd);
  TH1* h_triple = mix_histograms(h_double, h_edep, rnd);

  // Combine them with their expected rates
  h_edep->Scale(n_exp);
  TH1* h_combined = (TH1*) h_edep->Clone("combined");
  h_double->Scale(n_exp_double);
  h_triple->Scale(n_exp_triple);
  h_combined->Add(h_double);
  h_combined->Add(h_triple);

  // Plot the result
  h_edep->SetLineColor(kBlue);
  h_double->SetLineColor(kMagenta);
  h_triple->SetLineColor(kGreen);
  h_combined->SetLineColor(kRed);
  h_edep->SetLineWidth(2);
  h_double->SetLineWidth(2);
  h_triple->SetLineWidth(2);
  h_combined->SetLineWidth(2);
  h_edep->SetTitle(Form("Total Calo Energy;Energy (MeV);Rate / %.1g POT / %.1f MeV", npot, h_edep->GetBinWidth(1)));
  h_combined->Draw("hist");
  h_edep->Draw("hist same");
  h_double->Draw("hist same");
  h_triple->Draw("hist same");
  h_combined->Draw("hist same");
  TLegend legend(0.6, 0.7, 0.9, 0.9);
  legend.AddEntry(h_edep, "Single Edep", "l");
  legend.AddEntry(h_double, "Double Edep", "l");
  legend.AddEntry(h_triple, "Triple Edep", "l");
  legend.AddEntry(h_combined, "Combined", "l");
  legend.Draw();
  const double max_val = h_combined->GetMaximum();
  h_combined->GetYaxis()->SetRangeUser(0., 1.2*max_val);
  c.SaveAs(Form("%s/double_edep.png", fig_dir));
  h_combined->GetYaxis()->SetRangeUser(1e-10*max_val, 5.*max_val);
  c.SetLogy();
  c.SaveAs(Form("%s/double_edep_log.png", fig_dir));

  // Add the results to a histogram file
  TFile out_file(Form("%s/nts.user.double_edep.version.sequencer.root", run_dir), "RECREATE");
  h_edep->Write();
  h_double->Write();
  h_triple->Write();
  h_combined->Write();
  for(size_t index = 0; index < n_bkgs; ++index) {
    auto h = h_edeps[index];
    h->SetTitle(titles[index].Data());
    h->Write();
  }
  out_file.Close();

  return 0;
}
