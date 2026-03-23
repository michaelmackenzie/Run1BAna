// Test closure approximation weights

double closure_approx(const double energy, const double kmax) {
  if(energy < 0. || energy > kmax || kmax == 0.) return 0.;
  const double x = energy/kmax;
  const double p = 20.*(1. - 2.*x + 2.*x*x)*x*(1.-x)*(1.-x);
  return p;
}

void test_closure() {

  const double kmax = 90.1;
  const double emin =   0.;
  const double emax = 110.;

  const Long64_t nsamples = 1e6;
  TRandom3 rnd(90);
  TH1* h_gen = new TH1F("h_gen", "Generated spectrum", 200, 0., 100.);
  TH1* h_wt  = new TH1F("h_wt" , "Weighted spectrum" , 200, 0., 100.);

  for(Long64_t sample = 0; sample < nsamples; ++sample) {
    const double energy = rnd.Uniform(emin, emax);
    const double weight = (emax - emin)/kmax * closure_approx(energy, kmax);
    h_gen->Fill(energy);
    h_wt->Fill(energy, weight);
  }

  TCanvas* c = new TCanvas();
  gStyle->SetOptStat(1000011);
  h_gen->Draw();
  h_wt->Draw("sames hist");
  h_gen->GetYaxis()->SetRangeUser(0.1, 1.2*max(h_gen->GetMaximum(), h_wt->GetMaximum()));
}
