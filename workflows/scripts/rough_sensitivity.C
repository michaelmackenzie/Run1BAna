// Evaluate the sensitivity

//----------------------------------------------------------------------------------------------------------
// Main function
int rough_sensitivity(TString sig_file_name, double sig_eff,
                      TString bkg_file_name,
                      TString name,
                      const char* run_dir = ".") {


  // Open the histogram files
  TFile* sig_file = TFile::Open(sig_file_name, "READ");
  if(!sig_file || sig_file->IsZombie()) {
    std::cerr << __func__ <<  ": Error opening file: " << sig_file_name << std::endl;
    return 1;
  }
  TFile* bkg_file = TFile::Open(bkg_file_name, "READ");
  if(!bkg_file || bkg_file->IsZombie()) {
    std::cerr << __func__ <<  ": Error opening file: " << bkg_file_name << std::endl;
    return 1;
  }

  // N(POT) to assume
  const double npot = 7.2e17; // an example rate
  const double signal_br = (name == "flat_gamma") ? 1.e-6 : 1.e-9;

  // Get the total calo energy deposit per event histograms
  TH1* h_sig_edep = dynamic_cast<TH1*>(sig_file->Get("EDepAna/hist_0/total_calo_energy"));
  if(!h_sig_edep) {
    std::cerr << "Error: Signal histogram not found in file: " << sig_file_name << std::endl;
    return 1;
  }
  h_sig_edep->Scale(npot * signal_br * sig_eff / h_sig_edep->GetEntries());

  TH1* h_bkg_edep = dynamic_cast<TH1*>(bkg_file->Get("combined"));
  if(!h_bkg_edep) {
    std::cerr << "Error: Background histogram not found in file: " << bkg_file_name << std::endl;
    return 1;
  }
  h_bkg_edep->Scale(npot);

  // Plot the inputs
  const char* fig_dir = Form("%s/figures", run_dir);
  gSystem->Exec(Form("mkdir -p %s", fig_dir));
  gStyle->SetOptStat(0);
  TCanvas c("c", "Double Edep", 800, 600);

  // Plot the result
  h_sig_edep->SetLineColor(kBlue);
  h_bkg_edep->SetLineColor(kRed);
  h_sig_edep->SetLineWidth(2);
  h_bkg_edep->SetLineWidth(2);
  h_sig_edep->SetTitle(Form("Total Calo Energy;Energy (MeV);Rate / %.1g POT / %.1f MeV", npot, h_sig_edep->GetBinWidth(1)));

  h_sig_edep->Draw("hist");
  h_bkg_edep->Draw("hist same");

  TLegend legend(0.6, 0.7, 0.9, 0.9);
  legend.AddEntry(h_sig_edep, "Signal", "l");
  legend.AddEntry(h_bkg_edep, "Background", "l");
  legend.Draw();

  const double max_val = max(h_sig_edep->GetMaximum(), h_bkg_edep->GetMaximum());
  h_sig_edep->GetYaxis()->SetRangeUser(0., 1.2*max_val);
  c.SaveAs(Form("%s/%s_sig_vs_bkg.png", fig_dir, name.Data()));
  h_sig_edep->GetYaxis()->SetRangeUser(1e-12*max_val, 5.*max_val);
  c.SetLogy();
  c.SaveAs(Form("%s/%s_sig_vs_bkg.png", fig_dir, name.Data()));

  const int bin = h_sig_edep->GetMaximumBin();
  const double mpv = h_sig_edep->GetBinCenter(bin);
  const double height = h_sig_edep->GetBinContent(bin);
  int bin_1 = h_sig_edep->FindFirstBinAbove(height/2.);
  int bin_2 = h_sig_edep->FindLastBinAbove(height/2.);
  if(name == "flat_gamma") { // hard code this since the input isn't a delta line
    bin_1 = h_sig_edep->FindBin(75.);
    bin_2 = h_sig_edep->FindBin(90.);
  }

  const double x1 = h_sig_edep->GetBinCenter(bin_1);
  const double x2 = h_sig_edep->GetBinCenter(bin_2);
  const double fwhm = x2 - x1;

  const double bkg_rate = h_bkg_edep->Integral(bin_1, bin_2);
  const double sig_rate = h_sig_edep->Integral(bin_1, bin_2);
  const double sensitivity = (bkg_rate > 0.) ? sig_rate / std::sqrt(bkg_rate) : 0.;

  std::cout << "Signal MPV = " << mpv << " FWHM = " << fwhm
            << " signal rate = " << sig_rate
            << " background rate = " << bkg_rate
            << " s/sqrt(b) = " << sensitivity
            << std::endl;

  return 0;
}
