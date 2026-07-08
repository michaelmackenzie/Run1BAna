#include <TFile.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TEfficiency.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>

void plot_calo_trig_eff(const char* inputFilename = "nts.trig.ce_v40.root", TString tag = "ce_v40") {
  // 1. Setup styles and create target output directory
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  if(tag != "") tag = "_" + tag;
  TString outputDir = "figures/calo_trig_eff" + tag + "/";
  gSystem->mkdir(outputDir, true); // force creation of parent directories

  // 2. Open input ROOT file
  TFile* file = TFile::Open(inputFilename, "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error: Could not open file " << inputFilename << std::endl;
    return;
  }

  // 3. Fetch histograms (Adjust directory path strings if yours differ)
  TString basePath = "CaloTriggerEffAna/hist_1/"; // Offline cluster selection

  TH2F* h2_energy_vs_path = (TH2F*)file->Get(basePath + "energy_vs_path");
  TH1F* h_all_cluster_energy = (TH1F*)file->Get(basePath + "cluster_energy");

  if (!h2_energy_vs_path || !h_all_cluster_energy) {
    std::cerr << "Error: Found missing histograms inside " << basePath << std::endl;
    file->Close();
    return;
  }

  // 4. Create a canvas for drawing
  TCanvas* c1 = new TCanvas("c1", "Trigger Diagnostics", 800, 600);

  // Get number of bins on X axis (paths axis)
  int nBinsX = h2_energy_vs_path->GetNbinsX();

  // 5. Loop over every path bin inside the TH2F
  for (int binX = 1; binX <= nBinsX; ++binX) {
    // Retrieve the text label string set by your Fill(path.c_str(), ...) logic
    const char* pathName = h2_energy_vs_path->GetXaxis()->GetBinLabel(binX);

    // Skip default/empty bins that never held a string label
    if (!pathName || strlen(pathName) == 0) continue;

    TString sPathName(pathName);
    std::cout << "Processing plots for path: " << sPathName << std::endl;

    // Clean up path name string for filename safety (remove colons or slashes)
    TString safeName = sPathName;
    safeName.ReplaceAll(":", "_");
    safeName.ReplaceAll("/", "_");

    // ---- Plot 1: Absolute Deposited Energy Distribution ----
    c1->Clear();
    c1->SetLogy(0); // Set to 1 if you prefer logarithmic scale

    // Project a 1D slice of Y (Energy) for this specific path bin X
    TH1D* h_proj_energy = h2_energy_vs_path->ProjectionY(Form("proj_%s", safeName.Data()), binX, binX);
    h_proj_energy->SetTitle(Form("Cluster Energy passed by %s;Energy [MeV];Clusters", pathName));
    h_proj_energy->SetLineColor(kBlue);
    h_proj_energy->SetLineWidth(2);
    h_proj_energy->Draw("HIST");

    c1->SaveAs(outputDir + "energy_" + safeName + ".png");

    // ---- Plot 2: Efficiency Graph Relative to All Clusters ----
    c1->Clear();

    // Ensure the two histograms have identical bin alignment before passing to TEfficiency
    if (TEfficiency::CheckConsistency(*h_proj_energy, *h_all_cluster_energy)) {
      TEfficiency* eff = new TEfficiency(*h_proj_energy, *h_all_cluster_energy);
      eff->SetTitle(Form("Trigger Selection Efficiency: %s;Energy [MeV];Efficiency", pathName));
      eff->SetLineColor(kRed);
      eff->SetMarkerColor(kRed);
      eff->SetMarkerStyle(20);
      eff->SetMarkerSize(0.8);

      eff->Draw("AP"); // A: Draw axes, P: Draw markers with error bars

      // Adjust the Y-axis range strictly from 0 to 105% for standard viewing
      gPad->Update();
      if (eff->GetPaintedGraph()) {
        eff->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.0, 1.05);
      }
      gPad->Modified();

      c1->SaveAs(outputDir + "eff_" + safeName + ".png");
      delete eff;
    } else {
      std::cerr << "Warning: Bin mismatch between projection and total clusters for " << pathName << std::endl;
    }

    delete h_proj_energy;
  }

  // 6. Clean up memory allocations
  delete c1;
  file->Close();
  delete file;
  std::cout << "Done! All figures saved inside: " << outputDir << std::endl;
}
