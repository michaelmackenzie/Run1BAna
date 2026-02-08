// Compare signal momentum distribution in different selections
#include "Run1BAna/analysis/background_model.C"

int compare_p_dpst(TString process = "mumem", vector<int> selections = {20, 34, 35}) {
  // ensure selections are defined
  if(selections.empty()) return 0;

  const double ngen = get_dataset_info(process).ngen_;
  const double scale = (ngen > 0.) ? 1./ngen : 1.;

  // Retrieve the data
  vector<TH1*> hists;
  vector<double> integrals;
  double max_val = 0.;
  for(int selection : selections) {
    hists.push_back(get_background_hist(process, selection, process));
    TH1* h = hists.back();
    if(!h) return max(1, selection);

    // Scale it to per gen particle
    h->Scale(scale);
    const double sum = h->Integral();
    integrals.push_back(sum);

    // Normalize for plotting
    if(sum > 0.) h->Scale(1./sum);
    max_val = max(max_val, h->GetMaximum());

    // Set histogram style
    h->SetLineWidth(2);
    const int color = hists.size() + 1;
    h->SetLineColor(color);
  }

  // Plot the histograms
  gStyle->SetOptStat(0);
  TCanvas* c = new TCanvas("c","c", 1000, 700);
  TH1* haxis = hists[0];
  haxis->Draw("hist");
  for(size_t index = 1; index < hists.size(); ++index) hists[index]->Draw("same");

  haxis->SetTitle(""); haxis->SetXTitle("p (MeV/c)"); haxis->SetYTitle(Form("AU / %.2f MeV/c", haxis->GetBinWidth(1)));
  if(process == "mumem") haxis->GetXaxis()->SetRangeUser(100., 106.);
  if(max_val > 0.) haxis->GetYaxis()->SetRangeUser(1.e-4*max_val, 1.2*max_val);

  // Add a legend
  TLegend* leg = new TLegend(c->GetLeftMargin()+0.01, 1. - c->GetTopMargin() - 0.15,
                             1. - c->GetRightMargin()-0.01, 1. - c->GetTopMargin() - 0.01);
  leg->SetNColumns(min(3, (int) hists.size()));
  leg->SetFillColor(0); leg->SetFillStyle(0);
  leg->SetLineColor(0); leg->SetLineStyle(0);
  leg->SetTextSize(0.03);
  for(size_t index = 0; index < hists.size(); ++index) {
    leg->AddEntry(hists[index], Form("Selection %zu: #epsilon = %.3f", index, integrals[index]), "L");
  }
  leg->Draw();

  return 0;
}
