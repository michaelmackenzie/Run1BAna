// Plot the closure approximation

double closure_approx(const double energy, const double kmax) {
  if(energy < 0. || energy > kmax || kmax == 0.) return 0.;
  const double x = energy/kmax;
  const double p = 20.*(1. - 2.*x + 2.*x*x)*x*(1.-x)*(1.-x);
  return p;
}

void plot_closure() {

  const double kmax = 90.1;

  const int nbins = 1000;
  TRandom3 rnd(90);
  TH1* h_g = new TH1D("h_g", "Closure approximation;Energy (MeV);", nbins, 0., 100.);
  TH1* h_e = new TH1D("h_e", "Electron spectrum"    , nbins, 0., 100.);

  for(int bin = 1; bin < nbins; ++bin) {
    const double x = h_g->GetBinCenter(bin);
    const double p = closure_approx(x, kmax);
    const double w = h_g->GetBinWidth(bin);
    h_g->Fill(x,p*w);
    const int nsteps = 1000;
    for(int istep = 1; istep <= nsteps; ++istep) {
      h_e->Fill(x*istep/nsteps, p*w/nsteps);
    }
  }
  h_g->Scale(1./h_g->Integral()/h_g->GetBinWidth(1));
  h_e->Scale(1./h_e->Integral()/h_e->GetBinWidth(1));

  TCanvas* c = new TCanvas("c","c", 900, 600);
  c->SetRightMargin(0.03); c->SetLeftMargin(0.06);
  gStyle->SetOptStat(0);
  h_g->Draw("hist");
  h_e->Draw("hist sames");
  h_g->GetYaxis()->SetRangeUser(0., 1.1*max(h_g->GetMaximum(), h_e->GetMaximum()));
  h_g->SetLineWidth(3);
  h_e->SetLineWidth(3);
  h_g->SetLineColor(kRed);
  h_e->SetLineColor(kBlue);

  TLegend* leg = new TLegend(0.6, 0.75, 0.89, 0.89);
  leg->SetLineWidth(0);
  leg->AddEntry(h_g, "Photon spectrum");
  leg->AddEntry(h_e, "Electron spectrum");
  leg->Draw();

  const char* dir = "figures/closure_approx";
  gSystem->Exec(Form("mkdir -p %s", dir));
  c->SaveAs(Form("%s/closure.png", dir));
}
