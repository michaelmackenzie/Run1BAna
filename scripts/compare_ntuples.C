// Compare two ntuple version, normalizing per gen event

//-----------------------------------------------------------------------------------------
struct FileInfo {
    TFile* file;
    TString name;
    double nseen;
    double ngen;
    double scale;
};

//-----------------------------------------------------------------------------------------
void compare_hist(TString name, TString type, const int set, double xmin, double xmax, FileInfo& info_1, FileInfo& info_2, const TString& fig_dir) {
  TH1* h1 = (TH1*) info_1.file->Get(Form("Run1BAna/%s_%i/%s", type.Data(), set, name.Data()));
  TH1* h2 = (TH1*) info_2.file->Get(Form("Run1BAna/%s_%i/%s", type.Data(), set, name.Data()));
  if(!h1 || !h2) {
    cout << __func__ << ": Histogram " << name.Data() << "/" << type.Data() << "/" << set
         << " not found in one or both files\n";
    return;
  }
  h1->SetDirectory(0);
  h2->SetDirectory(0);
  h1->Scale(info_1.scale);
  h2->Scale(info_2.scale);

  if(xmin < xmax) {
    h1->GetXaxis()->SetRangeUser(xmin, xmax);
    h2->GetXaxis()->SetRangeUser(xmin, xmax);
  }
  const double max_val = max(h1->GetMaximum(), h2->GetMaximum());
  h1->GetYaxis()->SetRangeUser(0., 1.2*max_val);
  h1->GetYaxis()->SetTitle(Form("Events / Gen event / %.3g", h1->GetBinWidth(1)));

  TCanvas* c = new TCanvas(Form("c_%s_%s_%i", name.Data(), type.Data(), set), "", 1000, 900);
  TPad* pad1 = new TPad("pad1", "", 0., 0.3, 1., 1.); pad1->Draw();
  TPad* pad2 = new TPad("pad2", "", 0., 0., 1., 0.3); pad2->Draw();
  pad1->SetBottomMargin(0.02);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.3);

  pad1->cd();
  h1->SetLineColor(kRed);
  h2->SetLineColor(kBlue);
  h1->SetLineWidth(2);
  h2->SetLineWidth(2);
  h1->Draw("hist");
  h2->Draw("hist same");
  h1->GetXaxis()->SetLabelSize(0.);
  TLegend* leg = new TLegend(0.2, 0.8, 0.8, 0.89);
  leg->SetNColumns(2);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.06);
  leg->AddEntry(h1, info_1.name, "l");
  leg->AddEntry(h2, info_2.name, "l");
  leg->Draw();

  // Ratio plot
  pad2->cd();
  TH1* h_ratio = (TH1*) h1->Clone(Form("h_ratio_%s_%s_%i", name.Data(), type.Data(), set));
  h_ratio->Divide(h2);
  h_ratio->GetYaxis()->SetRangeUser(0.41, 1.59);
  h_ratio->GetYaxis()->SetTitle("Ratio");
  h_ratio->SetTitle("");
  h_ratio->GetXaxis()->SetTitle("");
  h_ratio->Draw("E1");
  h_ratio->GetXaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetTitleSize(0.1);
  h_ratio->GetYaxis()->SetTitleOffset(0.5);
  if(xmin < xmax) h_ratio->GetXaxis()->SetRangeUser(xmin, xmax);

  c->SaveAs(Form("%s/%s_%s_%i.png", fig_dir.Data(), name.Data(), type.Data(), set));
  h1->GetYaxis()->SetRangeUser(max_val/1.e5, 10.*max_val);
  pad1->SetLogy();
  c->SaveAs(Form("%s/%s_%s_%i_log.png", fig_dir.Data(), name.Data(), type.Data(), set));
  delete c;
}

//-----------------------------------------------------------------------------------------
int compare_ntuples(TString file_1, TString file_2,
                    TString name_1 = "v1", TString name_2 = "v2",
                    TString fig_tag = "comp") {
  // Open the two files
  TFile* f1 = TFile::Open(file_1, "READ");
  TFile* f2 = TFile::Open(file_2, "READ");
  if(!f1 || !f2) return 1;

  // Get the normalization info
  TH1* h_nseen_1 = (TH1*) f1->Get("Run1BAna/data/norm");
  TH1* h_nseen_2 = (TH1*) f2->Get("Run1BAna/data/norm");
  TH1* h_ngen_1  = (TH1*) f1->Get("Run1BAna/data/ngen");
  TH1* h_ngen_2  = (TH1*) f2->Get("Run1BAna/data/ngen");
  if(!h_nseen_1 || !h_ngen_1 || !h_nseen_2 || !h_ngen_2) {
    cout << __func__ << ": Normalization histograms not found in one or both files\n"
         << " File 1: nseen = " << (h_nseen_1 ? "found" : "not found") << ", ngen = " << (h_ngen_1 ? "found" : "not found") << endl
         << " File 2: nseen = " << (h_nseen_2 ? "found" : "not found") << ", ngen = " << (h_ngen_2 ? "found" : "not found") << endl;
    return 1;
  }
  const double nseen_1 = h_nseen_1->Integral();
  const double nseen_2 = h_nseen_2->Integral();
  const double ngen_1 = h_ngen_1->Integral();
  const double ngen_2 = h_ngen_2->Integral();
  const double scale_1 = (ngen_1 > 0.) ? 1. / ngen_1 : 0.;
  const double scale_2 = (ngen_2 > 0.) ? 1. / ngen_2 : 0.;

  cout << "File 1: nseen = " << nseen_1 << ", ngen = " << ngen_1 << ", scale = " << scale_1 << endl;
  cout << "File 2: nseen = " << nseen_2 << ", ngen = " << ngen_2 << ", scale = " << scale_2 << endl;

  FileInfo info_1 = {f1, name_1, nseen_1, ngen_1, scale_1};
  FileInfo info_2 = {f2, name_2, nseen_2, ngen_2, scale_2};

  // Make the figure directory
  TString fig_dir = "figures/compare_ntuples";
  if(fig_tag != "") fig_dir += "_" + fig_tag;
  gSystem->mkdir(fig_dir, true);
  gStyle->SetOptStat(0);

  // Compare histograms
  compare_hist("npot"       , "evt", 0, 0., 5.e7, info_1, info_2, fig_dir);
  compare_hist("ncombo_hits", "evt", 0, 1., -1. , info_1, info_2, fig_dir);
  compare_hist("ncalo_hits" , "evt", 0, 1., -1. , info_1, info_2, fig_dir);
  compare_hist("energy"     , "cls", 0, 0., 100., info_1, info_2, fig_dir);
  compare_hist("time"       , "cls", 0, 1., -1. , info_1, info_2, fig_dir);

  return 0;
}
