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
double nnt_sig_      = 0.;
double nnt_bkg_      = 0.;
double sig_skim_eff_ = 1.;
double norm_sig_     = 0.;
double norm_bkg_     = 0.;

//----- -------------------------------------------------------------------------
double getNSampled(TFile* f, int set = 0) {
  TH1* h_norm = (TH1*) f->Get(Form("Run1BAna/evt_%i/npot", set));
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
    return;
  }

  h_sig  = (TH1*) h_sig ->Clone(Form("%s_sig", name));
  h_bkg = (TH1*) h_bkg->Clone(Form("%s_bkg", name));
  const int norm_set = 60; // efficiencies relative to the base photon selection set
  const double eff_sig  = h_sig ->GetEntries() / getNSampled(f_sig, norm_set);
  const double eff_bkg = h_bkg->GetEntries() / getNSampled(f_bkg, norm_set);
  if(normalize) {
    if(eff_sig > 0.) h_sig ->Scale(1./normInRange(h_sig, x_min, x_max));
    if(eff_bkg > 0.) h_bkg->Scale(1./normInRange(h_bkg, x_min, x_max));
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

  h_sig ->SetLineColor(kBlue);
  h_bkg->SetLineColor(kRed);
  h_sig ->SetLineWidth(3);
  h_bkg->SetLineWidth(3);
  h_sig ->SetFillStyle(3004);
  h_bkg->SetFillStyle(3005);
  h_sig ->SetFillColor(kBlue);
  // h_bkg->SetFillColor(kRed);
  h_sig ->Draw("hist");
  h_bkg->Draw("hist same");

  const double max_sig  = h_sig ->GetMaximum();
  const double max_bkg = h_bkg->GetMaximum();
  const double max_val = std::max(max_sig, max_bkg);
  const double min_max = (max_sig <= 0.) ? max_bkg : (max_bkg <= 0.) ? max_sig : std::min(max_sig, max_bkg);
  h_sig->GetYaxis()->SetRangeUser(0., 1.2*max_val);

  if(min_max < 0.) {
    cout << "!!! " << name << "/" << type << "/" << set << ": Max(sig) = " << max_sig
         << " Max(bkg) = " << max_bkg << endl;
  }

  TLegend legend(0.6, 0.75, 0.9, 0.9);
  legend.AddEntry(h_sig , Form("Signal (eff = %.3g%%)"    , 100.*eff_sig));
  legend.AddEntry(h_bkg, Form("Background (eff = %.3g%%)", 100.*eff_bkg));
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.Draw();

  TString fig_name = Form("%s/%s_%s_%i%s", dir_.Data(), name, type, set, (normalize) ? "_norm" : "");
  c.SaveAs((fig_name + ".png").Data());

  double ymin = std::max(min_max*1.e-3, 1.e-6);
  double r = max_val/ymin;
  double factor = std::max(2., std::log10(r)*5.); // scale up proportional to orders of magnitude spanned
  double ymax = max_val*factor;
  h_sig->GetYaxis()->SetRangeUser(ymin, ymax);
  c.SetLogy();
  c.SaveAs((fig_name + "_log.png").Data());

  delete h_sig;
  delete h_bkg;

}

//-----------------------------------------------------------------------------------------------
void plot_gen_eff(TFile* f, int set) {
  TH1* h_gen = (TH1*) f->Get(Form("Run1BAna/sim_%i/energy_start", set)); // get generated energy distribution
  if(!h_gen) {
    std::cerr << "Cannot find energy_start histogram for set " << set << "\n";
    return;
  }
  h_gen = (TH1*) h_gen->Clone(Form("eff_%i", set));
  const int nsampled = ((TH1*) f->Get("Run1BAna/evt_0/npot"))->GetEntries();

  // Scale the generated energy histogram by the number of events and the provided scale factor
  h_gen->Scale(sig_skim_eff_ / nsampled);
  h_gen->Rebin(5); // Rebin to reduce fluctuations

  // Assume generation was flat between 50 and 110 MeV
  h_gen->Scale((110. - 50.) / h_gen->GetXaxis()->GetBinWidth(1));
  h_gen->SetTitle(Form("Generated energy distribution;Energy (MeV);Efficiency / %.1g MeV", h_gen->GetXaxis()->GetBinWidth(1)));

  TCanvas c("c","c", 1000, 800);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.05);

  h_gen->Draw("E1");
  h_gen->GetXaxis()->SetRangeUser(50., 110.);
  h_gen->SetLineWidth(5);
  h_gen->SetMarkerStyle(20);
  h_gen->SetMarkerSize(1.);
  h_gen->SetMarkerColor(kRed);
  h_gen->SetLineColor(kGray);

  c.SaveAs(Form("%s/photon_eff_%i.png", dir_.Data(), set));
}

//------------------------------------------------------------------------------
void plotRMCvsBkg(const char* filename_sig = "nts.mu2e.FlatGammaMix1BB-KL.Run1Bah_best_v1_4-001.root",
                  const char* filename_bkg = "nts.mu2e.NoPrimaryMix1BB-KL_skim_clusters.Run1Bah_best_v1_4-001.root",
                  // const char* filename_sig = "nts.mu2e.FlatGammaMixLowTriggerable-KL.Run1Baf_best_v1_4-000.root",
                  // const char* filename_bkg = "nts.mu2e.NoPrimaryMix1BB_skim_trig_clusters.Run1Bah_best_v1_4-000.root",
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
  // sig_skim_eff_ = (1263859. / 1949000000.); // digi dataset
  sig_skim_eff_ = (1257537. / 1940000000.); // mcs dataset
  const double norm_g  = (110. - 50.)/90.1 * sig_skim_eff_; // sample creation + filtering factors
  // const double norm_b  = 3880./100000.; // digi trigger cluster skim
  const double norm_b = 17551. / 1.e6; // mcs cluster skim
  norm_sig_ = nrmc*norm_g/nnt_sig_;
  norm_bkg_ = nevents*norm_b/nnt_bkg_;
  std::cout << "Norms: RMC = " << norm_sig_ << " Bkg = " << norm_bkg_ << std::endl;

  // Set up the figure directory and style
  dir_ = (tag) ? Form("figures/rmc_vs_bkg_%s", tag) : "figures/rmc_vs_bkg";
  gSystem->Exec(Form("mkdir -p %s", dir_.Data()));
  gStyle->SetOptStat(0);

  // Plot the histograms
  vector<int> sets = {0, 10, 60, 61, 62, 63, 67, 68, 69, 70};
  for(const int set : sets) {
    plot_gen_eff(f_sig, set);
    for(const bool normalize : {false, true}) {
      plot("energy", "cls", set, normalize, 1, 70., 120., f_sig, f_bkg);
      plot("time", "cls", set, normalize, 2, 400., 2000., f_sig, f_bkg);
      plot("radius", "cls", set, normalize, 1, 300., 700., f_sig, f_bkg);
      plot("disk", "cls", set, normalize, 1, 0., 2., f_sig, f_bkg);
      plot("frac_first_crystal", "cls", set, normalize, 1, 1., -1., f_sig, f_bkg);
      plot("frac_first_two_crystals", "cls", set, normalize, 1, 1., -1., f_sig, f_bkg);
      plot("second_moment", "cls", set, normalize, 5, 1., -1., f_sig, f_bkg);
      plot("t_var", "cls", set, normalize, 1, 0., 5., f_sig, f_bkg);
      plot("time_cluster_dt", "cls", set, normalize, 1, -50., 50., f_sig, f_bkg);
      plot("ncr", "cls", set, normalize, 1, 0., 10., f_sig, f_bkg);
      plot("nhits", "tcls", set, normalize, 4, 1., -1., f_sig, f_bkg);
      plot("nhits", "csms", set, normalize, 4, 1., -1., f_sig, f_bkg);
      plot("photon_id", "cls", set, normalize, 1, 0., 1., f_sig, f_bkg);
      // if(set == 60) {
      //   plot("A0", "csms", set, normalize, 4, 1., -1., f_sig, f_bkg);
      //   plot("A1", "csms", set, normalize, 4, 1., -1., f_sig, f_bkg);
      //   plot("B0", "csms", set, normalize, 4, 1., -1., f_sig, f_bkg);
      //   plot("B1", "csms", set, normalize, 4, 1., -1., f_sig, f_bkg);
      //   plot("d0", "csms", set, normalize, 4, 1., -1., f_sig, f_bkg);
      // }
    }
  }

}
