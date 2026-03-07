// Compare electron and muon distributions
#include "ConvAna/tools/utilities.C"
int Mode_ = 0;

int plot(TFile* f_1, TFile* f_2, TString hist, TString xtitle, TString unit,
         int set_1, int set_2, int rebin = 1, double xmin = 1., double xmax = -1., bool logy = false) {
  TH1* h_1 = (TH1*) f_1->Get(Form("Ana/ConvAna_ConvAna/Hist/trk_%i/%s", set_1, hist.Data()));
  TH1* h_2 = (TH1*) f_2->Get(Form("Ana/ConvAna_ConvAna/Hist/trk_%i/%s", set_2, hist.Data()));
  if(!h_1 || !h_2) return 1;
  if(h_1->Integral() <= 0. || h_2->Integral() <= 0.) return 1;

  h_1 = (TH1*) h_1->Clone("h_1");
  h_2 = (TH1*) h_2->Clone("h_2");

  h_1->Scale(1./h_1->Integral());
  h_2->Scale(1./h_2->Integral());

  if(rebin > 1) {
    h_1->Rebin(rebin);
    h_2->Rebin(rebin);
  }

  h_1->SetLineWidth(2); h_1->SetLineColor(kBlue); h_1->SetFillStyle(3004); h_1->SetFillColor(h_1->GetLineColor());
  h_2->SetLineWidth(2); h_2->SetLineColor(kRed ); h_2->SetFillStyle(3005); h_2->SetFillColor(h_2->GetLineColor());

  gStyle->SetOptStat(0);
  auto c = make_ratio_plot(Plot_t(h_1, {h_2}, xtitle, "", unit, xmin, xmax, 1., -1., 0., 2., logy, -1., -1.));
  if(!c) return 1;

  TLegend* leg = new TLegend(0.26, 0.85, 0.80, 0.92);
  leg->SetNColumns(2);
  leg->AddEntry(h_1, "Electron", "PLE");
  leg->AddEntry(h_2, "Muon"    , "FL");
  leg->SetLineWidth(0);
  c->cd();
  leg->Draw();

  c->SaveAs(Form("figures/m%i/%s_%i_%i.png", Mode_, hist.Data(), set_1, set_2));
  delete c;
  delete h_1;
  delete h_2;

  return 0;
}

/**
   Mode:
      0: CRY
      1: flat e- + flat mu-
 **/
int compare(int Mode = 0) {
  Mode_ = Mode;
  TFile* f_1 = nullptr;
  TFile* f_2 = nullptr;
  if(Mode == 0) {
    f_1 = TFile::Open("/exp/mu2e/data/users/mmackenz/conv_ana/histograms/ConvAna.cnv_ana.cry4ab1s5r0001.m1.hist", "READ");
    f_2 = TFile::Open("/exp/mu2e/data/users/mmackenz/conv_ana/histograms/ConvAna.cnv_ana.cry4ab1s5r0001.m1.hist", "READ");
  } else if(Mode == 1) {
    f_1 = TFile::Open("/exp/mu2e/data/users/mmackenz/conv_ana/histograms/ConvAna.cnv_ana.fele0b0s5r0001.m1.hist", "READ");
    f_2 = TFile::Open("/exp/mu2e/data/users/mmackenz/conv_ana/histograms/ConvAna.cnv_ana.fmum0b0s5r0001.m1.hist", "READ");
  }
  if(!f_1 || !f_2) return 1;

  gSystem->Exec(Form("[ ! -d figures/m%i ] && mkdir -p figures/m%i", Mode_, Mode_));
  int status(0);
  for(int region = 0; region <= 20; region += 10) {
    const int base_set = region + 150;
    status += plot(f_1, f_2, "tzslope"         , "dt/dz", "ns/mm"            , base_set, base_set+1,  2, -0.01, 0.03);
    status += plot(f_1, f_2, "tzslopesig"      , "dt/dz significance", ""    , base_set, base_set+1,  2,    -2,  10.);
    status += plot(f_1, f_2, "tzsloperatio"    , "dt/dz / expected", ""      , base_set, base_set+1,  0,   -1.,   3.);
    status += plot(f_1, f_2, "p"               , "q*p", "MeV/c"              , base_set, base_set+1,  5,  -150,  150);
    status += plot(f_1, f_2, "p_2"             , "p", "MeV/c"                , base_set, base_set+1, 40,   80,   110);
    status += plot(f_1, f_2, "cosTheta"        , "cos(#theta)", ""           , base_set, base_set+1,  5,   0.4,  0.9);
    status += plot(f_1, f_2, "fitCons_log"     , "log10(p(#chi^{2}))", ""    , base_set, base_set+1, 10,  -10.,   0.);
    status += plot(f_1, f_2, "fitMomErr"       , "#sigma(p)", "MeV/c"        , base_set, base_set+1,  2,    0.,  0.5);
    status += plot(f_1, f_2, "t0err"           , "#sigma(t)", "ns"           , base_set, base_set+1,  5,    0.,   5.);
    status += plot(f_1, f_2, "nActive"         , "N(active hits)", ""        , base_set, base_set+1,  2,   10.,  70.);
    status += plot(f_1, f_2, "nActiveFrac"     , "N(active hits)/N(hits)", "", base_set, base_set+1,  4,   0.5,   1.);
    status += plot(f_1, f_2, "ep"              , "E/p", ""                   , base_set, base_set+1,  1,    0.,  1.5);
    status += plot(f_1, f_2, "dt"              , "#deltat", "ns"             , base_set, base_set+1,  1,   -5.,  10.);
  }

  return status;
}
