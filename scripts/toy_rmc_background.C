// Plot the closure approximation
#include "Run1BAna/tools/functions.C"

//-------------------------------------------------------------------------------------
double closure_approx(const double energy, const double kmax) {
  if(energy < 0. || energy > kmax || kmax == 0.) return 0.;
  const double x = energy/kmax;
  const double p = 20.*(1. - 2.*x + 2.*x*x)*x*(1.-x)*(1.-x);
  return p;
}

//-------------------------------------------------------------------------------------
TH1* conversions_from_photons(TH1* h) {
  TH1* h_conv = (TH1*) h->Clone(Form("%s_conv", h->GetName()));
  h_conv->Reset();
  const int nbins = h_conv->GetNbinsX();
  const int nsteps = 1000; // conversion spectrum steps, assuming flat

  // Step through the photon spectrum
  for(int bin = 1; bin < nbins; ++bin) {
    const double x = h->GetBinCenter(bin);
    const double p = h->GetBinContent(bin);
    const double w = h->GetBinWidth(bin);
    for(int istep = 1; istep <= nsteps; ++istep) {
      h_conv->Fill(x*istep/nsteps, p*w/nsteps);
    }
  }
  h_conv->Scale(h->Integral() / h_conv->Integral()); // ensure rate is conserved
  return h_conv;
}

//-------------------------------------------------------------------------------------
void toy_rmc_background() {

  const char* dir = "figures/toy_rmc_background";
  gSystem->Exec(Form("mkdir -p %s", dir));

  const double kmax = 90.1;
  const double kend = 101.866;
  const double frac_kend = 0.01;

  const int nbins = 1020;
  TRandom3 rnd(90);
  TH1* h_1 = new TH1D("h_1", "RMC spectrum;Energy (MeV);", nbins, 0., 102.);
  TH1* h_2 = new TH1D("h_2", "RMC spectrum;Energy (MeV);", nbins, 0., 102.);

  for(int bin = 1; bin < nbins; ++bin) {
    const double x = h_1->GetBinCenter(bin);
    const double p_1 = closure_approx(x, kmax);
    const double p_2 = p_1*(1. - frac_kend) + frac_kend*closure_approx(x, kend);
    const double w = h_1->GetBinWidth(bin);
    h_1->Fill(x,p_1*w);
    h_2->Fill(x,p_2*w);
  }
  h_1->Scale(1./h_1->Integral()/h_1->GetBinWidth(1));
  h_2->Scale(1./h_2->Integral()/h_2->GetBinWidth(1));

  // Normalization info
  const double p_2_90 = h_2->Integral(h_2->FindBin(90.), nbins)*h_2->GetBinWidth(1);
  const double p_1_57 = h_1->Integral(h_1->FindBin(57.), nbins)*h_1->GetBinWidth(1);
  cout << "P(photon > 90 MeV) in two closures = " << p_2_90 << endl;

  const double nmuons = 3.e17;
  const double nrmc   = nmuons*0.609*(1.41e-5/p_1_57);
  const double p_conv = 1.e-3; // P(conversion) in Mu2e + reco

  // Set the norms and styles
  h_1->Scale(nrmc);
  h_2->Scale(nrmc);
  h_1->SetLineWidth(3);
  h_2->SetLineWidth(3);
  h_1->SetLineColor(kRed);
  h_2->SetLineColor(kBlue);

  // Create the conversion histograms
  TH1* h_1_e = conversions_from_photons(h_1); h_1_e->Scale(p_conv);
  TH1* h_2_e = conversions_from_photons(h_2); h_2_e->Scale(p_conv);

  TCanvas* c = new TCanvas("c","c", 900, 600);
  c->SetRightMargin(0.03); c->SetLeftMargin(0.06);
  gStyle->SetOptStat(0);

  TLegend* leg = new TLegend(0.6, 0.75, 0.89, 0.89);
  leg->SetLineWidth(0);
  leg->AddEntry(h_1, "Closure spectrum");
  leg->AddEntry(h_2, "Two closures spectrum");

  // Draw the photons
  h_1->Draw("hist");
  h_2->Draw("hist same");
  leg->Draw();
  h_1->SetTitle(Form("RMC photons per %.1e muon stops;Energy (MeV);", nmuons));

  double max_val = max(h_1->GetMaximum(), h_2->GetMaximum());
  h_1->GetYaxis()->SetRangeUser(0., 1.1*max_val);
  c->SaveAs(Form("%s/photon_spectrum.png", dir));
  h_1->GetYaxis()->SetRangeUser(1.e-6*max_val, 5.*max_val);
  h_1->GetXaxis()->SetRangeUser(60., 102.);
  c->SetLogy();
  c->SaveAs(Form("%s/photon_spectrum_log.png", dir));


  // Draw the conversions
  c->SetLogy(0);
  h_1_e->Draw("hist");
  h_2_e->Draw("hist same");
  leg->Draw();
  h_1_e->SetTitle(Form("RMC positrons per %.1e muon stops;Energy (MeV);", nmuons));

  // Example signal
  auto signal = landau_crystal_ball_tf1();
  signal->SetParameters(100., 92., 0.5, 10.4, 0.8, 0.4, 5., 10.);
  signal->SetLineColor(kGreen+1);
  signal->SetRange(85., 95.);
  signal->Draw("same");
  leg->AddEntry(signal, "#mu^{-}#rightarrowe^{+}", "L");

  max_val = max(h_1_e->GetMaximum(), h_2_e->GetMaximum());
  h_1_e->GetYaxis()->SetRangeUser(0., 1.1*max_val);
  c->SaveAs(Form("%s/electron_spectrum.png", dir));
  h_1_e->GetYaxis()->SetRangeUser(1.e-10*max_val, 5.*max_val);
  h_1_e->GetXaxis()->SetRangeUser(60., 102.);
  c->SetLogy();
  c->SaveAs(Form("%s/electron_spectrum_log.png", dir));
}
