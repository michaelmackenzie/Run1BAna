// Plot RMC vs. Bkg

#include "Run1BAna/analysis/physics.C"

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TError.h"

#include <algorithm>

TString dir_; // figure directory
double nnt_sig_ = 0.;
double nnt_bkg_ = 0.;
double norm_sig_ = 0.;
double norm_bkg_ = 0.;

//----- -------------------------------------------------------------------------
double getNSampled(TFile* f) {
  TH1* h_norm = (TH1*) f->Get("Run1BAna/evt_0/npot");
  if(!h_norm) {
    Error("getNSampled", "Could not retrieve normalization histogram!");
    return 0.;
  }
  const double entries = h_norm->GetEntries();
  if(entries <= 0.) {
    Error("getNSampled", "Normalization histogram has non-positive content!");
    return 0.;
  }
  return entries;
}

//------------------------------------------------------------------------------
double normInRange(TH1* h, double x_min, double x_max) {
  if(!h) return 0.;
  if(x_min >= x_max) return h->Integral(); // if invalid range, return total integral
  const int bin_min = max(1, min(h->GetXaxis()->GetNbins(), h->GetXaxis()->FindBin(x_min)));
  const int bin_max = max(1, min(h->GetXaxis()->GetNbins(), h->GetXaxis()->FindBin(x_max)));
  const double norm = h->Integral(bin_min, bin_max);
  return norm;
}

//------------------------------------------------------------------------------
void plot(const char* name, const char* type, const int set, const bool normalize,
          const int rebin, const double x_min, const double x_max,
          TFile* f_sig, TFile* f_bkg) {

  TH1* h_sig = (TH1*) f_sig->Get(Form("Run1BAna/%s_%i/%s", type, set, name));
  TH1* h_bkg = (TH1*) f_bkg->Get(Form("Run1BAna/%s_%i/%s", type, set, name));
  if(!h_sig || !h_bkg) {
    Error("plot", "Could not retrieve histograms! %s/%s/%i", name, type, set);
  }

  h_sig  = (TH1*) h_sig ->Clone(Form("%s_sig", name));
  h_bkg = (TH1*) h_bkg->Clone(Form("%s_bkg", name));
  const double eff_sig  = h_sig ->GetEntries() / nnt_sig_;
  const double eff_bkg = h_bkg->GetEntries() / nnt_bkg_;
  if(normalize) {
    h_sig ->Scale(1./normInRange(h_sig, x_min, x_max));
    h_bkg->Scale(1./normInRange(h_bkg, x_min, x_max));
  } else {
    h_sig ->Scale(norm_sig_ );
    h_bkg->Scale(norm_bkg_);
  }
  if(rebin > 1) {
    h_sig ->Rebin(rebin);
    h_bkg->Rebin(rebin);
  }
  if(x_min < x_max) {
    h_sig ->GetXaxis()->SetRangeUser(x_min, x_max);
    h_bkg->GetXaxis()->SetRangeUser(x_min, x_max);
  }
  TCanvas c("c","c", 1000, 800);
  c.SetLeftMargin(0.08);
  c.SetRightMargin(0.05);

  h_sig ->SetLineColor(kRed);
  h_bkg->SetLineColor(kBlue);
  h_sig ->SetLineWidth(3);
  h_bkg->SetLineWidth(3);
  h_sig ->SetFillStyle(3004);
  h_bkg->SetFillStyle(3005);
  h_sig ->SetFillColor(kRed);
  // h_bkg->SetFillColor(kBlue);
  h_sig ->Draw("hist");
  h_bkg->Draw("hist same");

  const double max_sig  = h_sig ->GetMaximum();
  const double max_bkg = h_bkg->GetMaximum();
  const double max_val = std::max(max_sig, max_bkg);
  const double min_max = std::min(max_sig, max_bkg);
  h_sig->GetYaxis()->SetRangeUser(0., 1.2*max_val);

  TLegend legend(0.6, 0.75, 0.9, 0.9);
  legend.AddEntry(h_sig , Form("Signal (eff = %.2g%%)"    , 100.*eff_sig));
  legend.AddEntry(h_bkg, Form("Background (eff = %.2g%%)", 100.*eff_bkg));
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.Draw();

  TString fig_name = Form("%s/%s_%s_%i%s", dir_.Data(), name, type, set, (normalize) ? "_norm" : "");
  c.SaveAs((fig_name + ".png").Data());

  h_sig->GetYaxis()->SetRangeUser(min_max*1.e-3, max_val*5.);
  c.SetLogy();
  c.SaveAs((fig_name + "_log.png").Data());

  delete h_sig;
  delete h_bkg;

}

//------------------------------------------------------------------------------
void plotRMCvsBkg(const char* filename_sig = "nts.mu2e.FlatGammaMixLowTriggerable-KL.Run1Baf_best_v1_4-000.root",
                  const char* filename_bkg = "nts.mu2e.NoPrimaryMix1BB_skim_trig_clusters.Run1Bah_best_v1_4-000.root",
                  const char* tag = "v03") {

  // Open the data files
  TFile* f_sig = TFile::Open(filename_sig, "READ");
  TFile* f_bkg = TFile::Open(filename_bkg, "READ");
  if (!f_sig || f_sig->IsZombie() ||
      !f_bkg || f_bkg->IsZombie()) {
    Error("plotRMCvsBkg", "Could not open files!");
    return;
  }

  // Get the normalization information
  nnt_sig_  = getNSampled(f_sig );
  nnt_bkg_ = getNSampled(f_bkg);
  if(nnt_sig_ <= 0. || nnt_bkg_ <= 0.) {
    Error("plotRMCvsBkg", "Invalid normalization!");
    return;
  }
  const double nevents = livetime_week_*duty_cycle_1bb_/1.695e-6; // N(events) in a week
  const double nmuons  = nevents*1.6e7*0.375*nmuons_per_pot_run1b_;
  const double nrmc    = nmuons*muon_capture_fraction_*br_rmc_/rmc_frac_57_; // N(RMC) assuming closure
  const double norm_c  = (102. - 70.)/90.1 * (1263859. / 1949000000.);
  norm_sig_ = nrmc*norm_c/nnt_sig_;
  norm_bkg_ = nevents/nnt_bkg_;

  // Set up the figure directory and style
  dir_ = (tag) ? Form("figures/rmc_vs_bkg_%s", tag) : "figures/rmc_vs_bkg";
  gSystem->Exec(Form("mkdir -p %s", dir_.Data()));
  gStyle->SetOptStat(0);

  // Plot the histograms
  vector<int> sets = {0, 10, 20, 25};
  for(const int set : sets) {
    for(const bool normalize : {false, true}) {
      plot("energy", "cls", set, normalize, 1, 75., 120., f_sig, f_bkg);
      plot("time", "cls", set, normalize, 2, 400., 2000., f_sig, f_bkg);
      plot("frac_first_crystal", "cls", set, normalize, 1, 1., -1., f_sig, f_bkg);
      plot("frac_first_two_crystals", "cls", set, normalize, 1, 1., -1., f_sig, f_bkg);
      plot("second_moment", "cls", set, normalize, 5, 1., -1., f_sig, f_bkg);
      plot("t_var", "cls", set, normalize, 1, 0., 5., f_sig, f_bkg);
      plot("ncr", "cls", set, normalize, 1, 0., 10., f_sig, f_bkg);
      plot("nhits", "tcls", set, normalize, 4, 1., -1., f_sig, f_bkg);
    }
  }

}
