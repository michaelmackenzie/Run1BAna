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

bool stack_bkgs_ = true;
double plot_livetime_ = 0.;
double plot_npot_ = 0.;
double plot_nmuons_ = 0.;
double br_sig_ = -1.;

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

//-------------------------------------------------------------------------------
TLatex* draw_info(const double scale = 0.75) {
  TLatex *logo = new TLatex();

  logo->SetNDC();
  const int ntens_pot = (plot_npot_ > 0.) ? int(std::log10(plot_npot_)) : 0;
  const float head_pot = plot_npot_ / std::pow(10.,std::max(0,ntens_pot));
  const int ntens_time = (plot_livetime_ > 0.) ? int(std::log10(plot_livetime_)) : 0;
  const float head_time = plot_livetime_ / std::pow(10.,std::max(0,ntens_time));
  TString lumistamp = Form("%.1f x 10^{%i} POT (1.5 kW); %.1f x 10^{%i} s",
                           head_pot, ntens_pot,
                           head_time, ntens_time);
  float textSize = 0.042 * 1.25 * scale;
  float extraOverTextSize  = 0.76;
  float extraTextSize = extraOverTextSize*textSize;

  const float x0(gPad->GetLeftMargin()+0.015), y0(1 - (gPad->GetTopMargin() - 0.017));
  logo->SetTextAlign(11);
  logo->SetTextSize(textSize);
  logo->SetTextFont(61);
  logo->DrawLatex(x0, y0, "Mu2e");
  logo->SetTextSize(0.042*scale);
  logo->SetTextFont(52);
  logo->DrawLatex(x0 + 0.08, y0,  "Simulation");
  logo->SetTextSize(extraTextSize);
  logo->SetTextFont(42);
  logo->SetTextAlign(31);
  if(br_sig_ > 0.) {
    const double rmue = br_sig_/muon_capture_fraction_;
    const int ntens_br = int(std::log10(rmue));
    const float head_br = rmue / std::pow(10.,ntens_br);
    lumistamp += Form("; R_{#mue} = %.1f x 10^{%i}", head_br, ntens_br);
  }
  if(plot_npot_ > 0.) logo->DrawLatex(1. - gPad->GetRightMargin(), y0, lumistamp);
  return logo;
}

//-------------------------------------------------------------------------------
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
double maxInRange(TH1* h, double x_min, double x_max) {
  if(!h) return 0.;
  const int nbins = h->GetNbinsX();
  const int bin_min = (x_min < x_max) ? max(1, min(nbins, h->GetXaxis()->FindBin(x_min))) : 1;
  const int bin_max = (x_min < x_max) ? max(1, min(nbins, h->GetXaxis()->FindBin(x_max))) : nbins;
  double val = 0.;
  for(int bin = bin_min; bin <= bin_max; ++bin) {
    const double binc = h->GetBinContent(bin);
    val = max(val, binc);
  }
  return val;
}

//------------------------------------------------------------------------------
TH1* significance_hist(TH1* h_sig, TH1* h_bkg, double x_min = 1., double x_max = -1.) {
  if(!h_sig || !h_bkg) return nullptr;
  TH1* h = (TH1*) h_sig->Clone(Form("%s_significance", h_sig->GetName()));
  h->Reset();
  h->SetFillStyle(0);
  const int nbins = h->GetNbinsX();
  const int bin_min = (x_min < x_max) ? max(1, min(nbins, h_sig->GetXaxis()->FindBin(x_min))) : 1;
  const int bin_max = (x_min < x_max) ? max(1, min(nbins, h_sig->GetXaxis()->FindBin(x_max))) : nbins;
  for(int bin = bin_min; bin <= bin_max; ++bin) {
    const double s = h_sig->GetBinContent(bin);
    const double b = h_bkg->GetBinContent(bin);
    const double sig = (b <= 0.) ? 0. : s/sqrt(b);
    h->SetBinContent(bin, sig);
    h->SetBinError  (bin, 0.);
  }
  return h;
}

//-------------------------------------------------------------------------------
TH1* smooth_tail(TH1* h, double x_min = 1., double x_max = -1., int rebin = 1) {
  gStyle->SetOptFit(0);
  TH1* h_s = (TH1*) h->Clone(Form("%s_smooth", h->GetName()));
  if(rebin > 1) h_s->Rebin(rebin);

  // First find where the signal falls off
  const int nbins = h->GetNbinsX();
  const int bin_min = (x_min < x_max) ? max(1, min(nbins, h_s->GetXaxis()->FindBin(x_min))) : 1;
  const int bin_max = (x_min < x_max) ? max(1, min(nbins, h_s->GetXaxis()->FindBin(x_max))) : nbins;
  const int bin_start = max(h_s->FindFirstBinAbove(0.), bin_min);

  int bin_last = bin_start;
  for(int bin = bin_start + 1; bin <= bin_max; ++bin) {
    if(h_s->GetBinContent(bin) <= 0.) break;
    bin_last = bin;
  }
  cout << bin_start << " " << bin_last << endl;

  if(bin_start > bin_max) return h_s;
  if(bin_last >= bin_max || bin_last < bin_min + 1) return h_s;

  const int bin_1 = max(bin_start, bin_last - 5);
  const int bin_2 = min(bin_max, bin_last + 3);
  const double x_1 = h_s->GetBinCenter(bin_1);
  const double x_2 = h_s->GetBinCenter(bin_2);
  TF1 f("f", "expo(0)", x_1, x_2);
  h_s->Fit(&f, "wRX0");

  for(int bin = bin_1; bin <= bin_max; ++bin) h_s->SetBinContent(bin, f.Eval(h_s->GetBinCenter(bin)));
  // TCanvas c; h->Draw("hist"); h_s->Draw("hist same"); h_s->SetLineColor(kRed); h->SetAxisRange(x_min, x_max, "X"); c.SetLogy(); c.SaveAs("tmp.png");
  return h_s;
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
  h_sig->SetTitle("");
  h_sig->Draw("hist");

  draw_info();
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
          const int rebin, const double x_min, const double x_max,
          const bool sig_plot = false, const bool smooth = false) {
  TH1* h_sig = nullptr;
  TH1* h_bkg = nullptr;
  TH1* h_bkg_no_calo_mu = nullptr;
  vector<TH1*> h_bkgs;
  THStack h_stack("h_stack", "Background stack");

  for(auto& process : processes_) {
    TH1* h_loc = (process.is_signal) ? h_sig : h_bkg;
    TH1* h = (TH1*) process.f->Get(Form("hist_%i/%s", set + process.set_offset, name));
    if(!h) {
      Error(__func__, "Could not retrieve histogram %s from process %s", name, process.name.Data());
      return;
    }
    if(smooth) h = smooth_tail(h, x_min, x_max, rebin);
    if(!h_loc) {
      h_loc = (TH1*) h->Clone(Form("h_%s_%i_%s", name, set, (process.is_signal) ? "sig" : "bkg"));
      h_loc->Scale(process.norm);
      if(process.is_signal) h_sig = h_loc;
      else                  h_bkg = h_loc;
    }
    else h_loc->Add(h, process.norm);

    // special case for background with no calo muons
    if(!process.is_signal && process.set_offset != 200) {
      if(!h_bkg_no_calo_mu) {
        h_bkg_no_calo_mu = (TH1*) h->Clone(Form("h_%s_%i_bkg_no_calo_mu", name, set));
        h_bkg_no_calo_mu->Scale(process.norm);
      } else h_bkg_no_calo_mu->Add(h, process.norm);
    }

    if(!process.is_signal) {
      h = (TH1*) h->Clone(Form("h_%s_%i_%s", name, set, process.name.Data()));
      h->Scale(process.norm);
      h->SetLineColor(process.color);
      h->SetTitle(process.name);
      h->SetLineWidth(1);
      if(stack_bkgs_) {
        h->SetFillColor(process.color);
        // h->SetLineColor(process.color+1);
        h->SetLineColor(kBlack);
      }
      h_bkgs.push_back(h);
    }
  }
  if(!h_sig || !h_bkg) {
    Error(__func__, "Could not retrieve histograms! %s/%i", name, set);
    return;
  }
  for(auto h : h_bkgs) h_stack.Add(h);

  if(normalize) {
    const double n_sig = normInRange(h_sig, x_min, x_max);
    const double n_bkg = normInRange(h_bkg, x_min, x_max);
    if(n_sig > 0.) h_sig->Scale(1./n_sig);
    if(n_bkg > 0.) {
      h_bkg->Scale(1./n_bkg);
      for(auto h : h_bkgs) h->Scale(1./n_bkg);
    }
  }
  if(!smooth && rebin > 1) {
    h_sig->Rebin(rebin);
    h_bkg->Rebin(rebin);
    h_bkg_no_calo_mu->Rebin(rebin);
    for(auto h : h_bkgs) h->Rebin(rebin);
  }
  if(x_min < x_max) {
    h_sig->GetXaxis()->SetRangeUser(x_min, x_max);
    h_bkg->GetXaxis()->SetRangeUser(x_min, x_max);
    for(auto h : h_bkgs) h->GetXaxis()->SetRangeUser(x_min, x_max);
  }

  TCanvas c("c","c", 1000, 800);
  TPad pad1("pad1", "pad1", 0., (sig_plot) ? 0.3 : 0.0, 1., 1.);
  TPad pad2("pad2", "pad2", 0., 0., 1., 0.3);
  pad1.SetLeftMargin(0.13);
  pad1.SetRightMargin(0.05);
  pad1.Draw();
  if(sig_plot) {
    pad2.SetLeftMargin(pad1.GetLeftMargin());
    pad2.SetRightMargin(pad1.GetRightMargin());
    pad1.SetBottomMargin(0.03);
    pad2.SetTopMargin(0.03);
    pad2.SetBottomMargin(0.35);
    pad2.Draw();
  }

  pad1.cd();
  TLegend legend(pad1.GetLeftMargin()+0.02, (sig_plot) ? 0.75 : 0.80, 1. - pad1.GetLeftMargin() - 0.02, 1. - pad1.GetTopMargin() - 0.01);
  legend.SetNColumns(3);
  legend.AddEntry(h_sig, "Signal", "F");
  if(!stack_bkgs_) legend.AddEntry(h_bkg, "Background", "F");
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
  h_sig->SetTitle("");
  h_sig->SetYTitle(Form("N(events) / %.1g", h_sig->GetBinWidth(1)));
  h_sig ->Draw("hist");
  if(stack_bkgs_) h_stack.Draw("hist same noclear");
  for(auto h : h_bkgs) {
    legend.AddEntry(h, h->GetTitle(), "F");
    if(!stack_bkgs_) h->Draw("hist same");
  }
  if(!stack_bkgs_) h_bkg->Draw("hist same");
  else             h_sig->Draw("hist same");


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

  // Make the significance plot if requested
  if(sig_plot) {
    pad2.cd();

    TH1* h_sig_full = significance_hist(h_sig, h_bkg, x_min, x_max); h_sig_full->SetName("significance_full");
    TH1* h_sig_cut  = significance_hist(h_sig, h_bkg_no_calo_mu, x_min, x_max);
    h_sig_cut->Draw("hist");
    h_sig_full->Draw("hist same");
    h_sig_full->SetLineColor(kRed);
    // h_sig_full->Draw("hist same");
    const double max_val = maxInRange(h_sig_cut, x_min, x_max);
    h_sig_cut->GetYaxis()->SetRangeUser(0., 1.2*max_val);
    if(x_min < x_max) h_sig_cut->GetXaxis()->SetRangeUser(x_min, x_max);

    const double text_size = 0.19;
    const double label_size = 0.13;
    const double y_offset = 0.35;
    h_sig->GetXaxis()->SetTitle("");
    h_sig->GetXaxis()->SetLabelSize(0.);
    h_sig->GetYaxis()->SetLabelSize(0.15*0.3/0.7);
    h_sig->GetYaxis()->SetTitleSize(text_size*0.3/0.7);
    h_sig->GetYaxis()->SetTitleOffset(y_offset*0.7/0.3);
    h_sig_cut->SetYTitle("S/#sqrt{B}");
    h_sig_cut->GetXaxis()->SetTitleSize(text_size);
    h_sig_cut->GetYaxis()->SetTitleSize(text_size);
    h_sig_cut->GetXaxis()->SetTitleOffset(0.70);
    h_sig_cut->GetYaxis()->SetTitleOffset(y_offset);
    h_sig_cut->GetXaxis()->SetLabelSize(label_size);
    h_sig_cut->GetYaxis()->SetLabelSize(label_size);

    TLegend* leg_2 = new TLegend(0.5, 0.9 - pad2.GetTopMargin(), 0.99 - pad2.GetRightMargin(), 0.99 - pad2.GetTopMargin());
    leg_2->SetNColumns(2); leg_2->SetLineWidth(0); leg_2->SetFillColor(0);
    leg_2->SetTextSize(0.10);
    leg_2->AddEntry(h_sig_full, "Full background", "L");
    leg_2->AddEntry(h_sig_cut, "No calo muons", "L");
    leg_2->Draw();

    pad1.cd();
  }

  draw_info((sig_plot) ? 1.1 : 0.75);
  TString fig_name = Form("%s/%s_%i%s", dir_.Data(), name, set, (normalize) ? "_norm" : "");
  c.SaveAs((fig_name + ".png").Data());

  double ymin = std::max(((normalize) ? 1.e-5 : min_max*1.e-3), 1.e-6);
  double r = max_val/ymin;
  double factor = std::max(2., std::log10(r)*20.); // scale up proportional to orders of magnitude spanned
  double ymax = max_val*factor;
  h_sig->GetYaxis()->SetRangeUser(ymin, ymax);
  pad1.SetLogy();
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

  h_sig->SetLineColor(kBlue);
  h_sig->SetLineWidth(3);
  h_bkg->SetLineWidth(3);
  h_sig->SetFillStyle(3004);
  h_bkg->SetFillStyle(kSolid);
  h_sig->SetFillColor(kBlue);
  h_bkg->SetLineColor(kBlack);
  h_bkg->SetFillColor(kRed);
  h_sig->Draw("hist");
  h_bkg->Draw("hist same");
  h_sig->Draw("hist same");

  const double max_sig  = h_sig ->GetMaximum();
  const double max_bkg = h_bkg->GetMaximum();
  const double max_val = std::max(max_sig, max_bkg);
  const double min_max = (max_sig <= 0.) ? max_bkg : (max_bkg <= 0.) ? max_sig : std::min(max_sig, max_bkg);
  h_sig->GetYaxis()->SetRangeUser(0., 1.2*max_val);
  h_sig->SetTitle("");

  if(min_max < 0.) {
    cout << "!!! " << name << "/" << set << ": Max(sig) = " << max_sig
         << " Max(bkg) = " << max_bkg << endl;
  }

  TLegend legend(0.6, 0.75, 0.9, 0.9);
  legend.AddEntry(h_sig, Form("Signal (eff = %.3g%%)"    , 100.*eff_sig), "F");
  legend.AddEntry(h_bkg, Form("Background (eff = %.3g%%)", 100.*eff_bkg), "F");
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.Draw();

  draw_info();
  TString fig_name = Form("%s/%s_%i%s", dir_.Data(), name, set, (normalize) ? "_norm" : "");
  c.SaveAs((fig_name + ".png").Data());

  double ymin = std::max(min_max*1.e-3, 1.e-6);
  double r = max_val/ymin;
  double factor = std::max(2., std::log10(r)*20.); // scale up proportional to orders of magnitude spanned
  double ymax = max_val*factor;
  h_sig->GetYaxis()->SetRangeUser(ymin, ymax);
  c.SetLogy();
  c.SaveAs((fig_name + "_log.png").Data());

  delete h_sig;
  delete h_bkg;

}

//-----------------------------------------------------------------------------------------------
void plot_gen_eff(TFile* f, int set) {
  TH1* h_gen = (TH1*) f->Get(Form("hist_%i/gen_energy_nowt", set)); // get generated energy distribution
  if(!h_gen) {
    std::cerr << "Cannot find energy_start histogram for set " << set << "\n";
    return;
  }
  h_gen = (TH1*) h_gen->Clone(Form("eff_%i", set));
  const int nsampled = getNSampled(f, -1);

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

  c.SaveAs(Form("%s/gen_eff_%i.png", dir_.Data(), set));
}
