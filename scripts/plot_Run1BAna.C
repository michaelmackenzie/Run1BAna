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

TH1* findHistBySubstring(TDirectory* dir, const std::string& substr) {
  if(!dir) return nullptr;
  TIter next(dir->GetListOfKeys());
  TKey* key;
  while((key = (TKey*)next())) {
    TString name = key->GetName();
    if(name.Contains(substr.c_str())) {
      TObject* obj = key->ReadObj();
      return dynamic_cast<TH1*>(obj);
    }
  }
  return nullptr;
}

void saveHistAsPNG(TH1* h,
                   const char* outdir, const char* fname,
                   int rebin = 1, double xmin = 1., double xmax = -1.,
                   const char* title=nullptr) {
  if(!h) {
    std::cerr << "Histogram " << fname << " not found.\n";
    return;
  }
  TCanvas c(fname, title ? title : fname, 800, 600);
  if(rebin > 1) h->Rebin(rebin);
  h->Draw("hist");
  h->SetLineColor(kBlue);
  h->SetLineWidth(3);
  h->SetFillColor(kBlue);
  h->SetFillStyle(3001);
  if(xmin < xmax) h->GetXaxis()->SetRangeUser(xmin, xmax);
  gPad->Update();
  TString out = TString::Format("%s/%s.png", outdir, fname);
  c.SaveAs(out.Data());
  c.SetLogy();
  out = TString::Format("%s/%s_log.png", outdir, fname);
  c.SaveAs(out.Data());
}

TH1* getHistogram(TFile* f, const char* dir, const char* name) {
  TH1* h = (TH1*) f->Get(Form("%s/%s", dir, name));
  if(h) h = (TH1*) h->Clone(Form("%s_%s", name, dir));
  return h;
}

void plot_Run1BAna(const char* infile = "run1b.root", const char* outdir = "plots") {
  gSystem->Exec(TString::Format("mkdir -p %s", outdir));

  TFile* f = TFile::Open(infile, "READ");
  if(!f || f->IsZombie()) {
    std::cerr << "Cannot open file: " << infile << "\n";
    return;
  }


  TH1* h_mc_digis  = getHistogram(f, "Run1BAna/evt_0", "mc_digis"  );
  TH1* h_combo     = getHistogram(f, "Run1BAna/evt_0", "combo_hits");
  TH1* h_calo      = getHistogram(f, "Run1BAna/evt_0", "calo_hits" );
  TH1* h_nclusters = getHistogram(f, "Run1BAna/evt_0", "nclusters" );
  TH1* h_hit_z     = getHistogram(f, "HitOriginAna/All events", "hit_z" );
  TH1* h_hit_t     = getHistogram(f, "HitOriginAna/All events", "hit_t" );
  TH1* h_hit_origin= getHistogram(f, "HitOriginAna/All events", "origin_type" );

  saveHistAsPNG(h_mc_digis,  outdir, "n_mc_digis"  , 5,     0., 5000., "N(MC digis)");
  saveHistAsPNG(h_combo,     outdir, "n_combo_hits", 5,     1.,   -1., "N(combo hits)");
  saveHistAsPNG(h_calo,      outdir, "n_calo_hits" , 1,     1.,   -1., "N(calo hits)");
  saveHistAsPNG(h_nclusters, outdir, "nclusters"   , 1,     1.,   -1., "N(calo clusters)");
  saveHistAsPNG(h_hit_z,     outdir, "hit_z"       , 1, -6000., 2000., "Straw digi z;z (mm)");
  saveHistAsPNG(h_hit_t,     outdir, "hit_t"       , 1,     1.,   -1., "Straw digi time;time (ns)");
  saveHistAsPNG(h_hit_origin,outdir, "hit_origin"  , 1,     1.,   -1., "Hit origin");

  // Cluster directory cls_0 -> energy
  for(int idir = 0; idir < 2; ++idir) {
    TDirectory* dc = (TDirectory*)f->Get(Form("Run1BAna/cls_%i", idir));
    if(!dc) {
      std::cerr << "Directory cls_0 not found in file.\n";
    }
    TH1* h_energy = findHistBySubstring(dc, "energy");
    saveHistAsPNG(h_energy, outdir, Form("cls_%i_energy", idir), 1, 0., 120., "Cluster energy");
    TH1* h_time = findHistBySubstring(dc, "time");
    saveHistAsPNG(h_time, outdir, Form("cls_%i_time", idir), 1, 300., 2000., "Cluster time");
  }

  f->Close();
  std::cout << "Plots saved to: " << outdir << "\n";
}
