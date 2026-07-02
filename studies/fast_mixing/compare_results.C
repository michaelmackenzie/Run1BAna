// Compare mixing results

//---------------------------------------------------------------------------
struct FileInfo_t {
  TFile* f = nullptr;
  Long64_t ngen = 1;
  TString version;
  int color = kBlue;
  FileInfo_t(TFile* f_, Long64_t ngen_, TString version_, int color_ = kBlue)
    : f(f_), ngen(ngen_), version(version_), color(color_) {}
};

//---------------------------------------------------------------------------
int plot(vector<FileInfo_t>& infos, TString hist, TString type, int set,
         double x_min = 1, double x_max = -1, double y_min = 1, double y_max = -1) {

  if(infos.size() < 2) {
    cerr << "Error: Not enough file infos provided" << endl;
    return 1;
  }

  // Retrieve the histograms
  vector<TH1*> hists;
  for(size_t index = 0; index < infos.size(); ++index) {
    FileInfo_t& info = infos[index];
    TH1* h = dynamic_cast<TH1*>(info.f->Get(Form("Run1BAna/%s_%i/%s", type.Data(), set, hist.Data())));
    if(!h) {
      cerr << "Error: Could not find histogram " << hist << " in file " << info.f->GetName() << endl;
      return 1;
    }
    h->SetDirectory(nullptr);
    h->SetName(Form("%s_%s_%i_%s", hist.Data(), type.Data(), set, info.version.Data()));
    h->Scale(1.0 / info.ngen);
    h->SetLineColor(info.color);
    h->SetLineWidth(2);
    hists.push_back(h);
  }

  // Draw the histograms
  TCanvas* c = new TCanvas(Form("c_%s_%s_%i", hist.Data(), type.Data(), set),
                           Form("c_%s_%s_%i", hist.Data(), type.Data(), set), 1000, 900);
  TPad* pad1 = new TPad("pad1", "pad1", 0., 0.3, 1., 1.); pad1->Draw();
  TPad* pad2 = new TPad("pad2", "pad2", 0., 0., 1., 0.3); pad2->Draw();
  pad1->SetBottomMargin(0.02);
  pad2->SetTopMargin(0.02);
  pad1->SetRightMargin(0.05);
  pad2->SetRightMargin(0.05);
  pad2->SetBottomMargin(0.3);

  pad1->cd();
  TLegend* legend = new TLegend(0.13, 0.75, 0.89, 0.89);
  legend->SetNColumns(min(4, (int)infos.size()));
  legend->SetBorderSize(0); legend->SetFillStyle(0); legend->SetLineWidth(0);
  legend->SetTextSize(0.04);
  TH1* haxis = nullptr;
  double max_value = hists[0]->GetMaximum();
  for(size_t index = 0; index < hists.size(); ++index) {
    TH1* h = hists[index];
    if(index == 0) {
      h->Draw("hist");
      haxis = h;
    } else {
      h->Draw("hist same");
      max_value = std::max(max_value, h->GetMaximum());
    }
    legend->AddEntry(h, infos[index].version, "L");
  }
  legend->Draw();
  haxis->GetXaxis()->SetTitleSize(0.);
  haxis->GetXaxis()->SetLabelSize(0.);
  haxis->GetYaxis()->SetTitle("Entries / N_{gen}");

  if(x_min > x_max) {
    x_min = haxis->GetXaxis()->GetXmin();
    x_max = haxis->GetXaxis()->GetXmax();
  }

  // Create a ratio plot
  pad2->cd();
  TH1* haxis_r = nullptr;
  vector<TH1*> hists_r;
  for(size_t index = 0; index < hists.size(); ++index) {
    TH1* h = hists[index];
    TH1* h_r = (TH1*)h->Clone(Form("%s_ratio", h->GetName()));
    h_r->Divide(hists[0]);
    h_r->SetLineColor(h->GetLineColor());
    h_r->SetLineWidth(h->GetLineWidth());
    if(index == 0) {
      h_r->Draw("hist");
      haxis_r = h_r;
    } else {
      h_r->Draw("hist same");
    }
    hists_r.push_back(h_r);
  }
  haxis_r->SetTitle("");
  haxis_r->GetXaxis()->SetTitleSize(0.1);
  haxis_r->GetXaxis()->SetTitleOffset(0.9);
  haxis_r->GetXaxis()->SetLabelSize(0.08);
  haxis_r->GetYaxis()->SetTitleSize(0.1);
  haxis_r->GetYaxis()->SetTitleOffset(0.5);
  haxis_r->GetYaxis()->SetLabelSize(0.08);
  haxis_r->GetYaxis()->SetNdivisions(505);
  haxis_r->GetYaxis()->SetRangeUser(0., 2.);
  haxis_r->GetYaxis()->SetTitle("Ratio");
  haxis_r->SetLineColor(0);

  pad1->cd();
  haxis  ->GetXaxis()->SetRangeUser(x_min, x_max);
  haxis_r->GetXaxis()->SetRangeUser(x_min, x_max);
  if(y_min < y_max) haxis->GetYaxis()->SetRangeUser(y_min, y_max);
  else              haxis->GetYaxis()->SetRangeUser(0, max_value * 1.2);
  c->SaveAs(Form("figures/compare_%s_%s_%i.png", hist.Data(), type.Data(), set));
  if(y_min >= y_max) haxis->GetYaxis()->SetRangeUser(max_value*1.e-5, max_value * 8);
  pad1->SetLogy();
  c->SaveAs(Form("figures/compare_%s_%s_%i_log.png", hist.Data(), type.Data(), set));

  delete c;
  delete legend;

  return 0;
}

//---------------------------------------------------------------------------
vector<FileInfo_t> getFileInfos(vector<TString> versions) {

  // Retrieve the data
  vector<FileInfo_t> infos;
  for(size_t index = 0; index < versions.size(); ++index) {
    TString version = versions[index];
    TString filename = Form("nts.mmackenz.fast_mixing.NoPrimary-reco.%s.root", version.Data());
    TFile* f = TFile::Open(filename);
    if(!f) {
      cerr << "Error: Could not open file " << filename << endl;
      return {};
    }
    TH1* h = dynamic_cast<TH1*>(f->Get("Run1BAna/data/ngen"));
    if(!h) {
      cerr << "Error: Could not find histogram in file " << filename << endl;
      return {};
    }
    infos.push_back(FileInfo_t(f, h->GetBinContent(1), version, index + 1 + (index > 3)));
  }
  return infos;
}

//---------------------------------------------------------------------------
int compare_results(vector<TString> versions = {"v0", "v1", "v2", "v3", "v6"}) {

  // Retrieve the data
  vector<FileInfo_t> infos = getFileInfos(versions);
  if(infos.size() != versions.size()) {
    cerr << "Error: Could not retrieve all file infos" << endl;
    return 1;
  }

  // Compare the data
  gSystem->Exec("mkdir -p figures");
  gStyle->SetOptStat(0);
  plot(infos, "energy", "cls", 1, 0., 120.);
  return 0;
}
