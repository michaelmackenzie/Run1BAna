// Plot the CAPHRI crystal energy spectrum

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TError.h"

#include <algorithm>

void plotCaphriEnergyOverlay(const char* filename,
                             bool normalize = false) {
  TFile* file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) {
    Error("plotCaphriEnergyOverlay", "Could not open file: %s", filename);
    return;
  }

  TH1* hAll    = dynamic_cast<TH1*>(file->Get("MuonGammaLineAna/info_0/hCaphriEnergy"));
  TH1* hSignal = dynamic_cast<TH1*>(file->Get("MuonGammaLineAna/info_1/hCaphriEnergy"));

  if (!hAll || !hSignal) {
    Error("plotCaphriEnergyOverlay",
          "Missing histogram(s). Expected info_0/hCaphriEnergy and info_1/hCaphriEnergy");
    file->Close();
    return;
  }

  TH1* hAllDraw = dynamic_cast<TH1*>(hAll->Clone("hCaphriEnergy_all_draw"));
  TH1* hSignalDraw = dynamic_cast<TH1*>(hSignal->Clone("hCaphriEnergy_signal_draw"));

  if (normalize) {
    if (hAllDraw->Integral() > 0.) hAllDraw->Scale(1.0 / hAllDraw->Integral());
    if (hSignalDraw->Integral() > 0.) hSignalDraw->Scale(1.0 / hSignalDraw->Integral());
    hAllDraw->GetYaxis()->SetTitle("Normalized entries");
  } else {
    hAllDraw->GetYaxis()->SetTitle("Events");
  }

  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("cCaphriEnergy", "CAPHRI hit energy overlay", 900, 650);

  hAllDraw->SetLineColor(kBlue + 1);
  hAllDraw->SetLineWidth(2);
  hAllDraw->SetTitle("CAPHRI Hit Energy;Energy (MeV);");

  hSignalDraw->SetLineColor(kRed + 1);
  hSignalDraw->SetLineWidth(2);


  hAllDraw->Draw("hist");
  hSignalDraw->Draw("hist same");

  TLegend* leg = new TLegend(0.60, 0.72, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hAllDraw, "All hits", "l");
  leg->AddEntry(hSignalDraw, "Signal-photon hits", "l");
  leg->Draw();

  double maxY = std::max(hAllDraw->GetMaximum(), hSignalDraw->GetMaximum());
  hAllDraw->GetXaxis()->SetRangeUser(0.5, 5.5);
  if(normalize) {
    hAllDraw->GetYaxis()->SetRangeUser(0., 1.15 * maxY);
    c->SaveAs("caphri_energy_norm.png");
    hAllDraw->GetYaxis()->SetRangeUser(0.001*maxY, 3. * maxY);
    c->SetLogy();
    c->SaveAs("caphri_energy_norm_log.png");
  } else {
    hAllDraw->GetYaxis()->SetRangeUser(0., 1.15 * maxY);
    c->SaveAs("caphri_energy.png");
    hAllDraw->GetYaxis()->SetRangeUser(min(0.01*hSignalDraw->GetMaximum(), 1.e-3*maxY), 5. * maxY);
    c->SetLogy();
    c->SaveAs("caphri_energy_log.png");
  }

  // file->Close();
}
