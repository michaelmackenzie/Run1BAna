// Plot Signal vs. Bkg

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

//-------------------------------------------------------------------------------
struct Process_t {
   TString name;
   TFile* f = nullptr;
   double norm = 1.;
   int set_offset = 0;
   bool is_signal = false;
   int color = kBlue;
};

//-------------------------------------------------------------------------------
vector<Process_t> processes_;

//----- -------------------------------------------------------------------------
double getIntegral(TFile* f, int set = 0) {
  if(set < 0) set = 0; // default -1 -> 0
  TH1* h_norm = (TH1*) f->Get(Form("hist_%i/cluster_energy", set));
  if(!h_norm) {
    Error(__func__, "Could not retrieve normalization histogram!");
    return 0.;
  }
  const double integral = h_norm->Integral(0, h_norm->GetNbinsX()+1);
  if(integral <= 0.) {
    Error(__func__, "Normalization histogram has non-positive integral!");
  }
  return integral;
}

//----- -------------------------------------------------------------------------
double getNSampled(TFile* f, int set = -1) {
  TH1* h_norm = (TH1*) ((set == -1) ? f->Get("norm") : f->Get(Form("hist_%i/cluster_energy", set)));
  if(!h_norm) {
    Error("getNSampled", "Could not retrieve normalization histogram for set %i!", set);
    return 0.;
  }
  const double entries = (set == -1) ? h_norm->Integral() : h_norm->GetEntries();
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
void plot_signal(TFile* f_sig, const char* name, const int set,
                 const int rebin, const double x_min, const double x_max) {
  TH1* h_sig = (TH1*) f_sig->Get(Form("hist_%i/%s", set, name));
  if(!h_sig) {
    Error(__func__, "Could not retrieve histogram %s from signal file!", name);
    return;
  }
  h_sig = (TH1*) h_sig->Clone(Form("h_%s_%i_signal", name, set));
  h_sig->Scale(norm_sig_);
  if(rebin > 1) h_sig->Rebin(rebin);
  if(x_min < x_max) h_sig->GetXaxis()->SetRangeUser(x_min, x_max);
  h_sig->SetLineColor(kBlue);
  h_sig->SetLineWidth(3);
  h_sig->SetFillStyle(3004);
  h_sig->SetFillColor(kBlue);
  TCanvas c("c","c", 1000, 800);
  c.SetLeftMargin(0.08);
  c.SetRightMargin(0.05);
  h_sig->Draw("hist");
  TString fig_name = Form("%s/%s_%i_signal", dir_.Data(), name, set);
  c.SaveAs((fig_name + ".png").Data());
  double max_val = h_sig->GetMaximum();
  h_sig->GetYaxis()->SetRangeUser(1.e-3*max_val, 1.2*max_val);
  c.SetLogy();
  c.SaveAs((fig_name + "_log.png").Data());
  delete h_sig;
}

//------------------------------------------------------------------------------
void plot(const char* name, const int set, const bool normalize,
          const int rebin, const double x_min, const double x_max) {
  TH1* h_sig = nullptr;
  TH1* h_bkg = nullptr;
  vector<TH1*> h_bkgs;

  for(auto& process : processes_) {
    TH1* h_loc = (process.is_signal) ? h_sig : h_bkg;
    TH1* h = (TH1*) process.f->Get(Form("hist_%i/%s", set + process.set_offset, name));
    if(!h) {
      Error(__func__, "Could not retrieve histogram %s from process %s", name, process.name.Data());
      return;
    }
    if(!h_loc) {
      h_loc = (TH1*) h->Clone(Form("h_%s_%i_%s", name, set, (process.is_signal) ? "sig" : "bkg"));
      h_loc->Scale(process.norm);
      if(process.is_signal) h_sig = h_loc;
      else                  h_bkg = h_loc;
    }
    else h_loc->Add(h, process.norm);
    if(!process.is_signal) {
      h = (TH1*) h->Clone(Form("h_%s_%i_%s", name, set, process.name.Data()));
      h->Scale(process.norm);
      h->SetLineColor(process.color);
      h->SetTitle(process.name);
      h_bkgs.push_back(h);
    }
  }
  if(!h_sig || !h_bkg) {
    Error(__func__, "Could not retrieve histograms! %s/%i", name, set);
    return;
  }

  if(normalize) {
    const double n_sig = normInRange(h_sig, x_min, x_max);
    const double n_bkg = normInRange(h_bkg, x_min, x_max);
    if(n_sig > 0.) h_sig->Scale(1./n_sig);
    if(n_bkg > 0.) {
      h_bkg->Scale(1./n_bkg);
      for(auto h : h_bkgs) h->Scale(1./n_bkg);
    }
  }
  if(rebin > 1) {
    h_sig->Rebin(rebin);
    h_bkg->Rebin(rebin);
    for(auto h : h_bkgs) h->Rebin(rebin);
  }
  if(x_min < x_max) {
    h_sig->GetXaxis()->SetRangeUser(x_min, x_max);
    h_bkg->GetXaxis()->SetRangeUser(x_min, x_max);
    for(auto h : h_bkgs) h->GetXaxis()->SetRangeUser(x_min, x_max);
  }

  TCanvas c("c","c", 1000, 800);
  c.SetLeftMargin(0.08);
  c.SetRightMargin(0.05);

  TLegend legend(0.5, 0.75, 0.9, 0.9);
  legend.SetNColumns(2);
  legend.AddEntry(h_sig, "Signal");
  legend.AddEntry(h_bkg, "Background");
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);

  h_sig->SetLineColor(kBlue);
  h_bkg->SetLineColor(kRed);
  h_sig->SetLineWidth(3);
  h_bkg->SetLineWidth(3);
  h_sig->SetFillStyle(3004);
  h_bkg->SetFillStyle(3005);
  h_bkg->SetLineStyle(kDashed);
  h_sig->SetFillColor(kBlue);
  // h_bkg->SetFillColor(kRed);
  h_sig ->Draw("hist");
  for(auto h : h_bkgs) {
    h->SetLineWidth(2);
    legend.AddEntry(h, h->GetTitle());
    h->Draw("hist same");
  }
  h_bkg->Draw("hist same");

  const double max_sig = h_sig->GetMaximum();
  const double max_bkg = h_bkg->GetMaximum();
  const double max_val = std::max(max_sig, max_bkg);
  const double min_max = (max_sig <= 0.) ? max_bkg : (max_bkg <= 0.) ? max_sig : std::min(max_sig, max_bkg);
  h_sig->GetYaxis()->SetRangeUser(0., 1.2*max_val);

  if(min_max < 0.) {
    cout << "!!! " << name << "/" << set << ": Max(sig) = " << max_sig
         << " Max(bkg) = " << max_bkg << endl;
  }

  legend.Draw();

  TString fig_name = Form("%s/%s_%i%s", dir_.Data(), name, set, (normalize) ? "_norm" : "");
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

//------------------------------------------------------------------------------
void plot(const char* name, const int set, const bool normalize,
          const int rebin, const double x_min, const double x_max,
          TFile* f_sig, TFile* f_bkg) {

  TH1* h_sig = (TH1*) f_sig->Get(Form("hist_%i/%s", set, name));
  TH1* h_bkg = (TH1*) f_bkg->Get(Form("hist_%i/%s", set, name));
  if(!h_sig || !h_bkg) {
    Error("plot", "Could not retrieve histograms! %s/%i", name, set);
    return;
  }

  h_sig  = (TH1*) h_sig ->Clone(Form("%s_sig", name));
  h_bkg = (TH1*) h_bkg->Clone(Form("%s_bkg", name));
  const int norm_set = 0; // efficiencies relative to the base photon selection set
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
    h_sig->Rebin(rebin);
    h_bkg->Rebin(rebin);
  }
  if(x_min < x_max) {
    h_sig->GetXaxis()->SetRangeUser(x_min, x_max);
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
    cout << "!!! " << name << "/" << set << ": Max(sig) = " << max_sig
         << " Max(bkg) = " << max_bkg << endl;
  }

  TLegend legend(0.6, 0.75, 0.9, 0.9);
  legend.AddEntry(h_sig , Form("Signal (eff = %.3g%%)"    , 100.*eff_sig));
  legend.AddEntry(h_bkg, Form("Background (eff = %.3g%%)", 100.*eff_bkg));
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.Draw();

  TString fig_name = Form("%s/%s_%i%s", dir_.Data(), name, set, (normalize) ? "_norm" : "");
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
