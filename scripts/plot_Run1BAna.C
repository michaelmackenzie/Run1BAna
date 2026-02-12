// ROOT macro to plot selected Run1BAna histograms
// Usage from shell: root -l -q 'plot_Run1BAna.C("file.root","outdir")'

#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TSystem.h>
#include <TString.h>
#include <iostream>
#include <string>


//-----------------------------------------------------------------------------------------------
TH1* getHistogram(TFile* f, const char* dir, const char* name) {
  TH1* h = (TH1*) f->Get(Form("%s/%s", dir, name));
  if(h) h = (TH1*) h->Clone(Form("%s_%s", name, dir));
  return h;
}

//-----------------------------------------------------------------------------------------------
void printPlot(TH1* h, const char* name,
               const char* outdir,
               int rebin = 1, double xmin = 1., double xmax = -1.,
               const char* title=nullptr) {
  TCanvas c(name, name, 800, 600);
  if(rebin > 1) h->Rebin(rebin);
  h->Draw("hist");
  h->SetLineColor(kBlue);
  h->SetLineWidth(3);
  h->SetFillColor(kBlue);
  h->SetFillStyle(3001);
  if(xmin < xmax) h->GetXaxis()->SetRangeUser(xmin, xmax);
  if(title) h->SetTitle(title);

  gPad->Update();
  c.SaveAs(Form("%s/%s.png", outdir, name));
  c.SetLogy();
  c.SaveAs(Form("%s/%s_log.png", outdir, name));
}

//-----------------------------------------------------------------------------------------------
void printPlot(TFile* f, const char* dir,
               const char* name, const char* outdir, const char* outname,
               int rebin = 1, double xmin = 1., double xmax = -1.,
               const char* title=nullptr) {
  TH1* h  = getHistogram(f, dir, name);
  if(!h) {
    std::cerr << "Histogram " << dir << "/" << name << " not found.\n";
    return;
  }
  printPlot(h, outname, outdir, rebin, xmin, xmax, title);
}

//-----------------------------------------------------------------------------------------------
void print_ratio_plot(TFile* f,
                      const char* dir_1, const char* name_1,
                      const char* dir_2, const char* name_2,
                      const char* outdir, const char* outname,
                      int rebin = 1, double xmin = 1., double xmax = -1.,
                      const char* title=nullptr) {
  TH1* h_1  = getHistogram(f, dir_1, name_1);
  TH1* h_2  = getHistogram(f, dir_2, name_2);
  if(!h_1) {
    std::cerr << "Histogram " << dir_1 << "/" << name_1 << " not found.\n";
    return;
  }
  if(!h_2) {
    std::cerr << "Histogram " << dir_2 << "/" << name_2 << " not found.\n";
    return;
  }
  TH1* h = (TH1*) h_1->Clone(Form("%s_ratio", h_1->GetName()));
  h->Divide(h_2);
  printPlot(h, outname, outdir, rebin, xmin, xmax, title);
}

//-----------------------------------------------------------------------------------------------
void plot_Run1BAna(const char* infile = "run1b.root", const char* outdir = "plots") {
  gSystem->Exec(TString::Format("mkdir -p %s", outdir));

  TFile* f = TFile::Open(infile, "READ");
  if(!f || f->IsZombie()) {
    std::cerr << "Cannot open file: " << infile << "\n";
    return;
  }

  printPlot(f, "Run1BAna/evt_0"         , "n_mc_digis"  , outdir,  "n_mc_digis"  , 5,     0., 5000., "N(MC digis)"              );
  printPlot(f, "Run1BAna/evt_0"         , "n_combo_hits", outdir,  "n_combo_hits", 5,     0., 5000., "N(combo hits)"            );
  printPlot(f, "Run1BAna/evt_0"         , "n_calo_hits" , outdir,  "n_calo_hits" , 1,     0.,  500., "N(calo hits)"             );
  printPlot(f, "Run1BAna/evt_0"         , "nclusters"   , outdir,  "nclusters"   , 1,     1.,   -1., "N(calo clusters)"         );
  printPlot(f, "Run1BAna/evt_0"         , "npot"        , outdir,  "npot"        , 1,     1.,   -1., "N(POT)"                   );
  printPlot(f, "HitOriginAna/All events", "hit_z"       , outdir,  "hit_z"       , 1, -6000., 2000., "Straw digi z;z (mm)"      );
  printPlot(f, "HitOriginAna/All events", "hit_t"       , outdir,  "hit_t"       , 1,     1.,   -1., "Straw digi time;time (ns)");
  printPlot(f, "HitOriginAna/All events", "origin_type" , outdir,  "origin_type" , 1,     1.,   -1., "Hit origin"               );

  // Cluster directory cls_0 -> energy
  for(int idir = 0; idir < 2; ++idir) {
    TString dir = Form("Run1BAna/cls_%i", idir);
    printPlot(f, dir.Data(), "energy", outdir, Form("energy_%i", idir), 1,   0.,  120., "Cluster energy;Energy (MeV)");
    printPlot(f, dir.Data(), "time"  , outdir, Form("time_%i"  , idir), 1, 300., 2000., "Cluster time;Time (ns)");
  }

  f->Close();
  std::cout << "Plots saved to: " << outdir << "\n";
}
