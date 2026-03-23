// Plot the CAPHRI crystal energy spectrum

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TError.h"

#include <algorithm>

TString dir_;


//------------------------------------------------------------------------------------------
TH1* fit_signal(TH1* h_sig) {
  TF1 f("signal_fit", "gaus");
  h_sig->Fit(&f);
  TH1* h = (TH1*) h_sig->Clone("h_signal_fit");
  h->Reset();
  for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
    const double xlow  = h_sig->GetBinLowEdge(bin);
    const double xhigh = h_sig->GetXaxis()->GetBinUpEdge(bin);
    if(xlow < 1.) continue;
    if(xhigh > 5.) break;
    h->SetBinContent(bin, f.Integral(xlow, xhigh));
  }
  return h;
}

//------------------------------------------------------------------------------------------
TH1* fit_background(TH1* h_bkg) {
  TF1 f("background_fit", "exp");
  h_bkg->Fit(&f);
  TH1* h = (TH1*) h_bkg->Clone("h_background_fit");
  h->Reset();
  for(int bin = 1; bin <= h->GetNbinsX(); ++bin) {
    const double xlow  = h_bkg->GetBinLowEdge(bin);
    const double xhigh = h_bkg->GetXaxis()->GetBinUpEdge(bin);
    h->SetBinContent(bin, f.Integral(xlow, xhigh));
  }
  return h;
}

//------------------------------------------------------------------------------------------
void caphriPlot(TFile* f, bool normalize) {
  TH1* hTot = (TH1*) f->Get("MuonGammaLineAna/info_0/hCaphriEnergy");
  TH1* hSig = (TH1*) f->Get("MuonGammaLineAna/info_2/hCaphriEnergy");

  if (!hTot || !hSig) {
    Error(__func__,
          "Missing histogram(s). Expected info_0/hCaphriEnergy and info_2/hCaphriEnergy");
    return;
  }

  hTot = (TH1*) hTot->Clone("h_all");
  hSig = (TH1*) hSig->Clone("h_signal");

  if (normalize) {
    if (hTot->Integral() > 0.) hTot->Scale(1.0 / hTot->Integral());
    if (hSig->Integral() > 0.) hSig->Scale(1.0 / hSig->Integral());
    hTot->GetYaxis()->SetTitle("Normalized entries");
  } else {
    hTot->GetYaxis()->SetTitle("Entries");
  }

  TCanvas* c = new TCanvas("cCaphriEnergy", "CAPHRI hit energy overlay", 900, 650);

  hTot->SetLineColor(kBlue + 1);
  hTot->SetLineWidth(2);
  hTot->SetTitle("CAPHRI Hit Energy;Energy (MeV);");

  hSig->SetLineColor(kRed + 1);
  hSig->SetLineWidth(2);


  hTot->Draw("hist");
  hSig->Draw("hist same");

  TLegend* leg = new TLegend(0.60, 0.72, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hTot, "All hits", "l");
  leg->AddEntry(hSig, "Signal-photon hits", "l");
  leg->Draw();

  double maxY = std::max(hTot->GetMaximum(), hSig->GetMaximum());
  hTot->GetXaxis()->SetRangeUser(0.5, 5.5);
  if(normalize) {
    hTot->GetYaxis()->SetRangeUser(0., 1.15 * maxY);
    c->SaveAs(Form("%s/caphri_energy_norm.png", dir_.Data()));
    hTot->GetYaxis()->SetRangeUser(0.001*maxY, 3. * maxY);
    c->SetLogy();
    c->SaveAs(Form("%s/caphri_energy_norm_log.png", dir_.Data()));
  } else {
    hTot->GetYaxis()->SetRangeUser(0., 1.15 * maxY);
    c->SaveAs(Form("%s/caphri_energy.png", dir_.Data()));
    hTot->GetYaxis()->SetRangeUser(min(0.01*hSig->GetMaximum(), 1.e-3*maxY), 5. * maxY);
    c->SetLogy();
    c->SaveAs(Form("%s/caphri_energy_log.png", dir_.Data()));
  }
  delete hTot;
  delete hSig;
}

//------------------------------------------------------------------------------------------
void fit_shapes(TFile* f) {
  TH1* hTot = (TH1*) f->Get("MuonGammaLineAna/info_0/hCaphriEnergy");
  TH1* hSig = (TH1*) f->Get("MuonGammaLineAna/info_2/hCaphriEnergy");

  if (!hTot || !hSig) {
    Error(__func__,
          "Missing histogram(s). Expected info_0/hCaphriEnergy and info_2/hCaphriEnergy");
    return;
  }

  const double nevents_per_second = 0.323/1.695e-6;
  const double nsampled = ((TH1*) f->Get("MuonGammaLineAna/info_0/hNCaloHits"))->GetEntries();
  const double norm = nevents_per_second / nsampled;

  hTot = (TH1*) hTot->Clone("h_all");
  hSig = (TH1*) hSig->Clone("h_signal");
  hTot->Scale(norm);
  hSig->Scale(norm);

  auto sig_fit = fit_signal(hSig);
  auto bkg_fit = fit_signal(hTot);
  cout << "N(signal photons / second) = " << sig_fit->Integral() << endl;

  TCanvas* c = new TCanvas("cCaphriEnergy", "CAPHRI hit energy overlay", 900, 650);

  bkg_fit->SetLineColor(kBlue + 1);
  bkg_fit->SetLineWidth(2);
  bkg_fit->SetTitle("CAPHRI Hit Energy per second;Energy (MeV);");

  sig_fit->SetLineColor(kRed + 1);
  sig_fit->SetLineWidth(2);


  bkg_fit->Draw("hist");
  sig_fit->Draw("hist same");

  TLegend* leg = new TLegend(0.60, 0.72, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(bkg_fit, "Background fit", "l");
  leg->AddEntry(sig_fit, "Signal-photon fit", "l");
  leg->Draw();

  double maxY = std::max(sig_fit->GetMaximum(), bkg_fit->GetMaximum());
  bkg_fit->GetXaxis()->SetRangeUser(0.5, 5.5);
  c->SaveAs(Form("%s/caphri_energy_fit.png", dir_.Data()));
  bkg_fit->GetYaxis()->SetRangeUser(min(0.01*hSig->GetMaximum(), 1.e-3*maxY), 5. * maxY);
  c->SetLogy();
  c->SaveAs(Form("%s/caphri_energy_fit_log.png", dir_.Data()));

  delete sig_fit;
  delete c;
}

//------------------------------------------------------------------------------------------
void observable(TFile* f, const char* name) {
  TH1* h = dynamic_cast<TH1*>(f->Get(Form("MuonGammaLineAna/info_0/%s", name)));

  if (!h) {
    cerr << __func__ << ": Histogram " << name << " not found\n";
    return;
  }

  h = dynamic_cast<TH1*>(h->Clone(Form("h_%s", name)));

  TCanvas c("c", "c", 900, 650);

  h->SetLineColor(kBlue + 1);
  h->SetLineWidth(2);

  h->Draw("hist");

  const double max_val = h->GetMaximum();
  h->GetYaxis()->SetRangeUser(0., 1.15*max_val);
  c.SaveAs(Form("%s/%s.png", dir_.Data(), name));
  h->GetYaxis()->SetRangeUser(0.001*max_val, 3. * max_val);
  c.SetLogy();
  c.SaveAs(Form("%s/%s_log.png", dir_.Data(), name));

  delete h;
}


//------------------------------------------------------------------------------------------
void plotCaphriEnergyOverlay(const char* filename = "nts.user.lumi_v06.version.sequencer.root",
                             bool normalize = false,
                             const char* tag = "v06") {
  TFile* file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) {
    Error("plotCaphriEnergyOverlay", "Could not open file: %s", filename);
    return;
  }

  dir_ = (tag) ? Form("figures/caphri_%s", tag) : "figures/caphri";
  gSystem->Exec(Form("mkdir -p %s", dir_.Data()));
  gStyle->SetOptStat(0);

  caphriPlot(file, normalize);

  fit_shapes(file);

  observable(file, "hNCaloHits"    );
  observable(file, "hCaloEnergy"   );
  observable(file, "hCaphriEnergy" );
  observable(file, "hNCaphriHits"  );
  observable(file, "hNTrackerHits" );
  observable(file, "hNTimeClusters");

  // file->Close();
}
