// Evaluate the sensitivity for Run 1A
// Assume signal = CE, background = DIO

//----------------------------------------------------------------------------------------------------------
// DIO theoretical spectrum
TH1* get_dio_spectrum() {
  TTree tree("t1","t1");
  TString table = "/exp/mu2e/app/users/mmackenz/run1b/Run1BAna/data/heeck_finer_binning_2016_szafron.tbl";
  int     nb    = 11000; //finer binning in this table
  double  bin   = 0.01;
  tree.ReadFile(table,"e/D:w/D");
  int n = tree.GetEntries();

  double emin = 0.;
  double emax = 110.;
  TH1* h_dio = new TH1D("h_dio","DIO spectrum",nb,emin,emax);

  double e, w;

  tree.SetBranchAddress("e",&e);
  tree.SetBranchAddress("w",&w);

  int prev_bin = 0;
  for (int i = 0; i < n; ++i) {
    tree.GetEntry(i);
    int ibin = h_dio->FindBin(e-bin/2.);
    if(prev_bin > 0 && prev_bin != ibin - 1) {
      cout << "Bin " << ibin << " entry " << i << " E = " << e
           << " but prev_bin = " << prev_bin << endl;
    }
    h_dio->SetBinContent(ibin,w);
    prev_bin = ibin;
  }
  h_dio->SetLineColor(kBlue);
  h_dio->SetLineWidth(2);
  h_dio->SetFillColor(0);
  h_dio->SetFillStyle(0);
  h_dio->Scale(1./(bin*h_dio->Integral()));
  return h_dio;
}

//----------------------------------------------------------------------------------------------------------
// Convolve a spectrum with a response
TH1* convolve(TH1* h_true, TH1* response) {
  const double eff = response->Integral()*response->GetBinWidth(1);
  TH1* h_reco = (TH1*) h_true->Clone(Form("%s_reco", h_true->GetName()));
  h_reco->Reset();
  for(int ibin = 1; ibin <= h_true->GetNbinsX(); ++ibin) {
    const double e_true = h_true->GetBinCenter (ibin);
    const double p_true = h_true->GetBinContent(ibin); // * h_true->GetBinWidth(ibin);
    for(int jbin = 1; jbin <= response->GetNbinsX(); ++jbin) {
      const double de = response->GetBinCenter(jbin);
      const double e_reco = e_true + de;
      const double p_resp = response->GetBinContent(jbin) * response->GetBinWidth(jbin);
      const double p = p_resp * p_true;
      h_reco->Fill(e_reco, p);
    }
  }
  h_reco->SetLineColor(kGreen-6);
  return h_reco;
}

//----------------------------------------------------------------------------------------------------------
// Get rough MPV and FWHM
void get_mpv_fwhm(TH1* h, double& mpv, double& fwhm) {
  int max_bin = h->GetMaximumBin();
  mpv = h->GetBinCenter(max_bin);
  const double max_val = h->GetBinContent(max_bin);
  const int bin_1 = h->FindFirstBinAbove(max_val/2.);
  const int bin_2 = max(bin_1, h->FindLastBinAbove(max_val/2.));
  const double x1 = h->GetBinLowEdge(bin_1);
  const double x2 = h->GetXaxis()->GetBinUpEdge(bin_2);
  fwhm = x2 - x1;
}


//----------------------------------------------------------------------------------------------------------
// Tracker resolution
TH1* trk_resolution() {
  const int nbins = 5000;
  TH1* res = new TH1D("response", "response", nbins, -5, 5);
  const double mu = 0.;
  const double sigma = 0.2;
  for(int bin = 1; bin <= nbins; ++bin) {
    const double p = ROOT::Math::gaussian_pdf(res->GetBinCenter(bin), sigma, mu);
    res->SetBinContent(bin, p);
  }
  return res;
}

//----------------------------------------------------------------------------------------------------------
// Main function
int rough_run1a_sensitivity(TString sig_file_name, double sig_eff,
                            const char* run_dir = ".") {


  // Open the histogram files
  TFile* sig_file = TFile::Open(sig_file_name, "READ");
  if(!sig_file || sig_file->IsZombie()) {
    std::cerr << __func__ <<  ": Error opening file: " << sig_file_name << std::endl;
    return 1;
  }

  // Setup output
  const char* fig_dir = Form("%s/figures", run_dir);
  gSystem->Exec(Form("mkdir -p %s", fig_dir));
  gStyle->SetOptStat(0);

  // N(POT) to assume
  const double npot = 1.e18; // an example rate
  const double signal_br = 1.e-13/0.609; // R_mue = 1.e-9
  const double mean_pot = 1.6e7; // 1BB
  const double nevents = npot / mean_pot;
  const double seconds = nevents * 1.695e-6; // seconds of On-Spill
  const double cosmic_rate_second = 2e4 / 1.1e7; // rough rate per second per MeV/c
  const double cosmic_rate = cosmic_rate_second * seconds; // rate per MeV/c

  // Get the signal energy and energy loss distribution before the tracker, for events with 10 MeV in the calo
  TH1* h_sig    = (TH1*) sig_file->Get("EDepAna/hist_2/trk_front_energy");
  TH1* response = (TH1*) sig_file->Get("EDepAna/hist_2/trk_front_energy_diff");
  if(!response || !h_sig) {
    std::cerr << "Error: Signal trk energy/ediff histograms not found in file: " << sig_file_name << std::endl;
    return 1;
  }
  h_sig->Scale(npot * signal_br * sig_eff / h_sig->GetEntries());
  response->Scale(sig_eff / response->GetEntries() / response->GetBinWidth(1));
  TH1* res = trk_resolution();
  h_sig = convolve(h_sig, res);
  h_sig->SetLineWidth(2);
  h_sig->SetLineColor(kBlue);
  double mpv, fwhm;
  get_mpv_fwhm(h_sig, mpv, fwhm);

  TCanvas c;
  response->Draw("hist");
  c.SaveAs(Form("%s/response.png", fig_dir));

  res->Draw("hist");
  c.SaveAs(Form("%s/res.png", fig_dir));

  // Get the DIO contribution
  TH1* dio = get_dio_spectrum();
  dio->Scale(0.39*sig_eff*npot); // DIO rate
  // TH1* dio_resp = convolve(dio, res);
  // TH1* dio_resp = convolve(dio, response);
  TH1* dio_resp = convolve(convolve(dio, response), res);
  dio->Draw("hist");
  dio_resp->Draw("hist same");
  c.SetLogy();
  c.SaveAs(Form("%s/dio.png", fig_dir));
  dio_resp->Rebin(h_sig->GetBinWidth(1) / dio->GetBinWidth(1));
  dio_resp->SetLineColor(kRed+1);

  // Get the Cosmic contribution
  TH1* cosmic = (TH1*) h_sig->Clone("cosmic");
  cosmic->Reset();
  for(int bin = 1; bin <= cosmic->GetNbinsX(); ++bin) {
    cosmic->SetBinContent(bin, cosmic_rate * cosmic->GetBinWidth(bin));
  }
  cosmic->SetLineColor(kGreen - 6);

  // Plot the signal + background
  h_sig->Draw("hist");
  dio_resp->Draw("hist same");
  cosmic->Draw("hist same");
  h_sig->SetTitle(Form("Signal vs. Background;Energy (MeV);Rate / %.1f MeV", h_sig->GetBinWidth(1)));
  h_sig->GetXaxis()->SetRangeUser(min(95., mpv - 1.5*fwhm), max(105., mpv + 1.5*fwhm));
  h_sig->GetYaxis()->SetRangeUser(1.e-3, 1.e4);

  TLegend leg(0.15, 0.8, 0.85, 0.89);
  leg.SetNColumns(3);
  leg.SetTextSize(0.05);
  leg.SetLineWidth(0);
  leg.SetFillColor(0);
  leg.AddEntry(h_sig, "Signal");
  leg.AddEntry(dio_resp, "DIO");
  leg.AddEntry(cosmic, "Cosmic");
  leg.Draw();

  c.SaveAs(Form("%s/sig_vs_bkg.png", fig_dir));


  // Evaluate the rough sensitivity
  double x_1(0.), x_2(0.), signal_rate(0.), dio_bkg(0.), cosmic_bkg(0.), bkg_rate(0.), sensitivity = -1.;

  for(int ibin = 1; ibin <= h_sig->GetNbinsX(); ++ibin) {
    for(int jbin = ibin; jbin <= h_sig->GetNbinsX(); ++jbin) {
      const double x_1_l = h_sig->GetBinCenter(ibin);
      const double x_2_l = h_sig->GetBinCenter(jbin);
      // if(x_1_l < mpv - 1.5*fwhm) continue;
      if(x_2_l < mpv) continue;
      const double signal_rate_l = h_sig   ->Integral(h_sig   ->FindBin(x_1_l), h_sig   ->FindBin(x_2_l));
      const double dio_bkg_l     = dio_resp->Integral(dio_resp->FindBin(x_1_l), dio_resp->FindBin(x_2_l));
      const double cosmic_bkg_l  = cosmic  ->Integral(cosmic  ->FindBin(x_1_l), cosmic  ->FindBin(x_2_l));
      const double bkg_rate_l    = dio_bkg_l + cosmic_bkg_l;
      if(signal_rate_l <= 0. || bkg_rate_l <= 0.) continue;
      const double sensitivity_l = signal_rate_l / std::sqrt(bkg_rate_l);
      if(sensitivity_l > sensitivity) {
        x_1 = x_1_l;
        x_2 = x_2_l;
        signal_rate = signal_rate_l;
        dio_bkg     = dio_bkg_l;
        cosmic_bkg  = cosmic_bkg_l;
        bkg_rate    = bkg_rate_l;
        sensitivity = sensitivity_l;
      }
    }
  }

  printf("Signal box = [%.1f, %.1f] MeV/c, signal = %.2g, dio = %.2g, cosmic = %.2g --> bkg = %.2g, S/sqrt(B) = %.3g\n",
         x_1, x_2, signal_rate, dio_bkg, cosmic_bkg, bkg_rate, sensitivity);
  return 0;
}
