// Make an approximation of the double DIO background

//---------------------------------------------------------------------------
// DIO theoretical spectrum
TH1* get_dio_spectrum() {
  TTree tree("t1","t1");
  TString table = "Run1BAna/data/heeck_finer_binning_2016_szafron.tbl";
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

//---------------------------------------------------------------------------
// Convolve the spectrum with calo response
TH1* double_dio(TH1* h_true) {
  // TH1* h_double = (TH1*) h_true->Clone(Form("%s_double", h_true->GetName()));
  TH1* h_double = new TH1D(Form("%s_double", h_true->GetName()), "Double DIO", 1100, 0., 110.); // use coarser binning
  const double e_max = h_true->GetXaxis()->GetBinUpEdge(h_true->GetNbinsX());
  h_double->Reset();
  for(int ibin = 1; ibin <= h_true->GetNbinsX(); ++ibin) {
    const double e_1 = h_true->GetBinCenter (ibin);
    const double p_1 = h_true->GetBinContent(ibin) * h_true->GetBinWidth(ibin);
    for(int jbin = 1; jbin <= h_true->GetNbinsX(); ++jbin) {
      const double e_2 = h_true->GetBinCenter (jbin);
      const double p_2 = h_true->GetBinContent(jbin) * h_true->GetBinWidth(jbin);
      h_double->Fill(e_1 + e_2, p_1*p_2);
    }
  }
  h_double->Scale(1./(h_double->Integral()*h_double->GetBinWidth(1)));
  h_double->SetLineColor(kAtlantic);
  return h_double;
}

//---------------------------------------------------------------------------
// Calo response: 5% resolution, ignore peak offset
double response(double e_true, double e_reco) {
  if(e_true <= 0. || e_reco < 0.) return 0.;

  const double u = (e_reco - e_true) / e_true;

  // Just do 5% Gaussian for now, no offset
  const double sigma = 0.05*e_true;
  const double mean  = 0.;
  return ROOT::Math::gaussian_pdf((e_reco - e_true), sigma, mean);
}

//---------------------------------------------------------------------------
// Convolve the spectrum with calo response
TH1* convolve(TH1* h_true) {
  const double de = std::min(0.1, h_true->GetBinWidth(1)); // resolution width step
  const double e_max = h_true->GetXaxis()->GetBinUpEdge(h_true->GetNbinsX());
  TH1* h_reco = (TH1*) h_true->Clone(Form("%s_reco", h_true->GetName()));
  h_reco->Reset();
  for(int ibin = 1; ibin <= h_true->GetNbinsX(); ++ibin) {
    const double e_true = h_true->GetBinCenter (ibin);
    const double p_true = h_true->GetBinContent(ibin) * h_true->GetBinWidth(ibin);
    double e_reco = de/2.;
    while(e_reco < e_max) {
      const double p_reco = response(e_true, e_reco);
      h_reco->Fill(e_reco, p_reco*p_true);
      e_reco += de;
    }
  }
  h_reco->SetLineColor(kGreen-6);
  return h_reco;
}

//---------------------------------------------------------------------------
void double_dio_mc(const double mean_muons_per_event = 20000., const double nmuons = 1e16) {

  // Get the PDFs
  TH1* h_dio = get_dio_spectrum();
  TH1* h_dio_double = double_dio(h_dio);

  // Rebin the DIO to match the double DIO
  h_dio->Rebin(10); h_dio->Scale(1./10.);

  // Get rate information
  const double r_dio           = 0.39; // DIO rate per stop
  const double dz_to_calo      = 4500.; // assuming 4.5 m from disk 0
  const double geom_overlap    = pow(34.,2)/(4.*M_PI*pow(dz_to_calo, 2)); // within a 34x34 mm area at 4.5 m away
  const double geom_acceptance =  (M_PI*(pow(675.,2) - pow(325.,2))) / (4.*M_PI*pow(dz_to_calo, 2)); // rough area of disk in acceptance
  const double time_overlap    = 0.0372; // roughly probability of two DIOs being coincident from dio_time_overlap.C
  const double acc_overlap     = geom_overlap*time_overlap; // P(overlap) given an accepted DIO (assuming N(muons) << 1/p_accept)
  const double p_overlap       = acc_overlap*r_dio*(mean_muons_per_event-1.); // P(overlap) given an accepted DIO (assuming N(muons) << 1/p_accept)
  h_dio->Scale(geom_acceptance*(1. - p_overlap));
  h_dio_double->Scale(geom_acceptance*p_overlap);

  // Print out the rate info
  cout << "-----------------------------------------------------------------\n"
       << "P(acceptance)      = " << geom_acceptance << endl
       << "P(overlap space)   = " << geom_overlap << endl
       << "P(overlap time)    = " << time_overlap << endl
       << "P(overlap)         = " << acc_overlap << endl
       << "P(any DIO overlap) = " << p_overlap << endl;

  // Convolve the histograms
  TH1* h_dio_reco        = convolve(h_dio);
  TH1* h_dio_double_reco = convolve(h_dio_double);
  h_dio_double_reco->SetLineColor(kOrange);

  // Print out the integrals
  cout << "-----------------------------------------------------------------\n"
       << "DIO                = " << h_dio            ->Integral()*h_dio            ->GetBinWidth(1) << endl
       << "DIO reco           = " << h_dio_reco       ->Integral()*h_dio_reco       ->GetBinWidth(1) << endl
       << "Double DIO         = " << h_dio_double     ->Integral()*h_dio_double     ->GetBinWidth(1) << endl
       << "Double DIO reco    = " << h_dio_double_reco->Integral()*h_dio_double_reco->GetBinWidth(1) << endl
       << "-----------------------------------------------------------------\n";

  // Scale the contributions to the given normalization
  h_dio            ->Scale(nmuons*r_dio);
  h_dio_reco       ->Scale(nmuons*r_dio);
  h_dio_double     ->Scale(nmuons*r_dio);
  h_dio_double_reco->Scale(nmuons*r_dio);

  // Rebin the distributions
  h_dio            ->Rebin(10);
  h_dio_reco       ->Rebin(10);
  h_dio_double     ->Rebin(10);
  h_dio_double_reco->Rebin(10);

  // Make the total reco histogram
  TH1* h_reco = (TH1*) h_dio_reco->Clone("h_reco");
  h_reco->Add(h_dio_double_reco);
  h_reco->SetLineColor(kRed);
  h_reco->SetLineStyle(kDashed);

  // Draw the histograms
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c", "c", 1200, 700);
  h_dio->Draw("hist");
  h_dio_reco->Draw("hist same");
  h_dio_double->Draw("hist same");
  h_dio_double_reco->Draw("hist same");
  h_reco->Draw("hist same");

  // Configure the axes
  const double ymax = h_dio->GetMaximum();
  h_dio->GetYaxis()->SetRangeUser(1.e-6, 100.*ymax);
  c->SetLogy();

  // Plot titles
  h_dio->SetTitle("");
  h_dio->SetXTitle("Energy (MeV)");
  h_dio->SetYTitle(Form("Events / %.2g MeV", h_dio->GetBinWidth(1)));

  // Make a legend
  TLegend* leg = new TLegend(0.59, 0.67, 0.89, 0.89);
  leg->SetFillColor(0); leg->SetFillStyle(0); leg->SetLineWidth(0); leg->SetLineColor(0); // don't show the box
  leg->SetNColumns(2);
  leg->SetTextSize(0.032);
  leg->AddEntry(h_dio            , "DIO");
  leg->AddEntry(h_dio_reco       , "DIO reco");
  leg->AddEntry(h_dio_double     , "Double DIO");
  leg->AddEntry(h_dio_double_reco, "Double DIO reco");
  leg->AddEntry(h_reco           , "Total reco");
  leg->Draw();

  // Make labels
  TLatex *logo = new TLatex();
  logo->SetNDC();
  const float x0(gPad->GetLeftMargin()+0.015), y0(1 - (gPad->GetTopMargin() - 0.017));
  logo->SetTextSize(0.032);
  logo->SetTextFont(42);
  logo->SetTextAlign(31);
  logo->DrawLatex(1. - gPad->GetRightMargin(), y0, Form("N(muon stops) = %.1e", nmuons));
}
