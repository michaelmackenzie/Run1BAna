#ifndef __RUN1BANA_TOOLS_UTILITIES__
#define __RUN1BANA_TOOLS_UTILITIES__

//------------------------------------------------------------------------------------------------------------
// General plot info
struct Plot_t {
  Plot_t(TH1* h,
         TString xtitle = "", TString ytitle = "", TString unit = "",
         double xmin = 1., double xmax = -1.,
         double ymin = 1., double ymax = -1.,
         double rmin = 0.5, double rmax = 1.5, bool logy = false,
         double npot = -1., double livetime = -1.) : h_(h), hists_({}),
                                                     xtitle_(xtitle), ytitle_(ytitle), unit_(unit),
                                                     xmin_(xmin), xmax_(xmax), ymin_(ymin), ymax_(ymax),
                                                     rmin_(rmin), rmax_(rmax), logy_(logy),
                                                     npot_(npot), livetime_(livetime)
  {}

  Plot_t(vector<TH1*> hists,
         TString xtitle = "", TString ytitle = "", TString unit = "",
         double xmin = 1., double xmax = -1.,
         double ymin = 1., double ymax = -1.,
         double rmin = 0.5, double rmax = 1.5, bool logy = false,
         double npot = -1., double livetime = -1.) : h_(nullptr), hists_(hists),
                                                     xtitle_(xtitle), ytitle_(ytitle), unit_(unit),
                                                     xmin_(xmin), xmax_(xmax), ymin_(ymin), ymax_(ymax),
                                                     rmin_(rmin), rmax_(rmax), logy_(logy),
                                                     npot_(npot), livetime_(livetime)
  {}

  Plot_t(TH1* h, vector<TH1*> hists,
         TString xtitle = "", TString ytitle = "", TString unit = "",
         double xmin = 1., double xmax = -1.,
         double ymin = 1., double ymax = -1.,
         double rmin = 0.5, double rmax = 1.5, bool logy = false,
         double npot = -1., double livetime = -1.) : h_(h), hists_(hists),
                                                     xtitle_(xtitle), ytitle_(ytitle), unit_(unit),
                                                     xmin_(xmin), xmax_(xmax), ymin_(ymin), ymax_(ymax),
                                                     rmin_(rmin), rmax_(rmax), logy_(logy),
                                                     npot_(npot), livetime_(livetime)
  {}

  TH1* h_;
  vector<TH1*> hists_;
  TString xtitle_;
  TString ytitle_;
  TString unit_;
  double xmin_;
  double xmax_;
  double ymin_;
  double ymax_;
  double rmin_;
  double rmax_;
  bool   logy_;
  double npot_;
  double livetime_;
};

//------------------------------------------------------------------------------------------------------------
TLatex* Mu2e_lumi(const bool is_data, const double npot = -1., const double livetime = -1., double scale = 1.) {
  TLatex *logo = new TLatex();

  logo->SetNDC();
  const int ntens_pot = (npot > 0.) ? int(std::log10(npot)) : 0;
  const float head_pot = npot / std::pow(10.,std::max(0,ntens_pot));
  const int ntens_time = (livetime > 0.) ? int(std::log10(livetime)) : 0;
  const float head_time = livetime / std::pow(10.,std::max(0,ntens_time));
  TString lumistamp = Form("%.1f x 10^{%i} POT (1BB); %.1f x 10^{%i} s",
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
  logo->DrawLatex(x0 + 0.08, y0,  (is_data) ? "Preliminary" : "Simulation");
  logo->SetTextSize(extraTextSize);
  logo->SetTextFont(42);
  logo->SetTextAlign(31);
  if(npot > 0.) logo->DrawLatex(1. - gPad->GetRightMargin(), y0, lumistamp);
  return logo;
}

//------------------------------------------------------------------------------------------------------------
//static function to delete all objects in a canvas
static Int_t Empty_Canvas(TCanvas* c) {
  if(!c) return 0;
  TList* list = c->GetListOfPrimitives();
  if(!list) return 1;
  for(auto o : *list) {
    if(o->InheritsFrom("TPad")) {
      auto pad = (TPad*) o;
      if(!pad) continue;
      TList* pad_list = pad->GetListOfPrimitives();
      for(auto h : *pad_list) {
        if(h->InheritsFrom("THStack")) {
          THStack* hstack = (THStack*) h;
          TList* hist_list = hstack->GetHists();
          for(auto hl : *hist_list) {
            delete hl;
          }
        }
        delete h;
      }
    }
    delete o;
  }
  if(c) delete c;
  return 0;
}

//------------------------------------------------------------------------------------------------------------
TH1* trim_hist(TH1* h, double xmin, double xmax) {
  if(xmin > xmax) return (TH1*) h->Clone(Form("%s_trim", h->GetName()));
  const int bin_low(h->FindBin(xmin+1.e-6)), bin_high(h->FindBin(xmax-1.e-6));
  TH1* h_trim = new TH1F(Form("%s_trim", h->GetName()), h->GetTitle(), (bin_high - bin_low +1), xmin, xmax);
  for(int ibin = bin_low; ibin <= bin_high; ++ibin) {
    const int bin = ibin - bin_low + 1;
    h_trim->SetBinContent(bin, h->GetBinContent(ibin));
    h_trim->SetBinError  (bin, h->GetBinError  (ibin));
  }
  return h_trim;
}

//------------------------------------------------------------------------------------------------------------
double min_in_range(TH1* h, double xmin, double xmax, bool use_errors = false, double cut_off = 0.) {
  const int bin_low((xmin > xmax) ? 1 : h->FindBin(xmin+1.e-6)), bin_high((xmin > xmax) ? h->GetNbinsX() : h->FindBin(xmax-1.e-6));
  double min_val = max(h->GetMaximum(), cut_off);
  for(int ibin = bin_low; ibin <= bin_high; ++ibin) {
    double binc = h->GetBinContent(ibin);
    if(binc <= cut_off) continue;
    if(use_errors) binc -= h->GetBinError(ibin);
    min_val = min(min_val, binc);
  }
  return min_val;
}

//------------------------------------------------------------------------------------------------------------
double max_in_range(TH1* h, double xmin, double xmax, bool use_errors = false) {
  const int bin_low((xmin > xmax) ? 1 : h->FindBin(xmin+1.e-6)), bin_high((xmin > xmax) ? h->GetNbinsX() : h->FindBin(xmax-1.e-6));
  double max_val = h->GetMinimum();
  for(int ibin = bin_low; ibin <= bin_high; ++ibin) {
    double binc = h->GetBinContent(ibin);
    if(use_errors) binc += h->GetBinError(ibin);
    max_val = max(max_val, binc);
  }
  return max_val;
}

//------------------------------------------------------------------------------------------------------------
TGraphAsymmErrors* errors_from_th1(TH1* hnom, TH1* hup, TH1* hdown = nullptr) {
  if(!hnom || !hup) {
    cout << __func__ << ": Inputs aren't defined!\n";
    return nullptr;
  }

  const int nbins(hnom->GetNbinsX());
  if(hup->GetNbinsX() != nbins || (hdown && hdown->GetNbinsX() != nbins)) {
    cout << __func__ << ": Shifted histograms have different binning!\n";
    return nullptr;
  }

  double xs[nbins], ys[nbins], xerrs[nbins], yups[nbins], ydowns[nbins];
  for(int ibin = 1; ibin <= nbins; ++ibin) {
    const int index = ibin - 1;
    const double x = hnom->GetBinCenter(ibin);
    const double y = hnom->GetBinContent(ibin);
    const double width = hnom->GetBinWidth(ibin);
    const double up = hup->GetBinContent(ibin) - y;
    const double down = (hdown) ? hdown->GetBinContent(ibin) - y : 0.;
    xs    [index] = x;
    ys    [index] = y;
    xerrs [index] = width / 2.;
    yups  [index] = max(0., max(up, down));
    ydowns[index] = max(0., -1.*min(up, down));
  }
  auto g = new TGraphAsymmErrors(nbins, xs, ys, xerrs, xerrs, ydowns, yups);
  g->SetMarkerStyle(hnom->GetMarkerStyle());
  g->SetMarkerSize(hnom->GetMarkerSize());
  g->SetMarkerColor(hnom->GetMarkerColor());
  g->SetFillStyle(3001);
  g->SetFillColor(kRed);
  return g;
}

//------------------------------------------------------------------------------------------------------------
// FIXME: Not finished
TCanvas* make_ratio_plot(Plot_t plot) {
  TH1* h = plot.h_; auto list = plot.hists_;
  if(!h || list.empty()) return nullptr;

  // Retrieve the plot info
  double xmin(plot.xmin_), xmax(plot.xmax_), ymin(plot.ymin_), ymax(plot.ymax_), rmin(plot.rmin_), rmax(plot.rmax_);
  TString xtitle(plot.xtitle_), ytitle(plot.ytitle_), unit(plot.unit_);
  bool logy = plot.logy_;

  // Configure the canvas
  TGaxis::SetExponentOffset(-0.06, 0.01, "Y"); // ensure the exponent doesn't hit label info
  TCanvas* c = new TCanvas("c", "c", 1200, 1000);
  TPad* pad1 = new TPad("pad1", "pad1", 0., 0.3, 1., 1.0); pad1->Draw();
  TPad* pad2 = new TPad("pad2", "pad2", 0., 0.0, 1., 0.3); pad2->Draw();
  pad1->SetBottomMargin(0.03); pad1->SetTopMargin(0.10); pad1->SetLeftMargin(0.12); pad1->SetRightMargin(0.05);
  pad2->SetBottomMargin(0.35); pad2->SetTopMargin(0.03); pad2->SetLeftMargin(pad1->GetLeftMargin()); pad2->SetRightMargin(pad1->GetRightMargin());
  pad1->SetFillColor(0); pad1->SetTickx(0); pad1->SetTicky(0);
  pad2->SetFillColor(0); pad2->SetTickx(0); pad2->SetTicky(0);

  // Check the axis limits
  if(xmax < xmin) {
    xmin = h->GetXaxis()->GetXmin();
    xmax = h->GetXaxis()->GetXmax();
  }
  xmin = std::max(xmin, h->GetXaxis()->GetXmin());
  xmax = std::min(xmax, h->GetXaxis()->GetXmax());

  if(ymax < ymin) {
    ymax = max_in_range(h, xmin, xmax, true);
    for(auto o : list) ymax = max(ymax, max_in_range(o, xmin, xmax));
    if(logy) {
      ymax = 50.*ymax;
      ymin = ymax/1.e5;
    } else {
      ymax = 1.2*ymax;
      ymin = 0.;
    }
  }

  // Draw the histograms
  pad1->cd();
  h->Draw("E1");
  for(auto o : list) o->Draw("hist same");
  h->Draw("E1 same"); // bring to the front

  // Configure the labels/fonts
  h->SetTitle("");
  h->GetYaxis()->SetTitle((ytitle == "") ? Form("Events / (%.2g%s)", h->GetXaxis()->GetBinWidth(1), (unit == "") ? "" : (" " + unit).Data()) : ytitle);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelOffset(0.01);
  h->GetYaxis()->SetTitleSize(0.06);
  h->GetYaxis()->SetTitleOffset(0.95);
  h->GetYaxis()->SetMaxDigits(3);
  h->GetXaxis()->SetLabelSize(0);
  h->GetXaxis()->SetLabelOffset(0.007);
  h->GetXaxis()->SetTitleSize(0.07);
  h->GetXaxis()->SetTitleOffset(0.80);
  h->SetTitle("");
  h->SetXTitle("");
  h->GetXaxis()->SetRangeUser(xmin, xmax);
  h->GetYaxis()->SetRangeUser(ymin, ymax);

  pad2->cd();

  // get ratios
  vector<TH1*> ratios;
  TH1* haxis = nullptr;
  for(auto num : list) {
    TH1* hr = (TH1*) num->Clone(Form("%s_r", num->GetName()));
    hr->Divide(h);
    if(ratios.empty()) {hr->Draw("E1"); haxis = hr;}
    else                hr->Draw("E1 same");
    ratios.push_back(hr);
  }

  if(haxis) {
    haxis->GetYaxis()->SetRangeUser(rmin, rmax);
    haxis->GetXaxis()->SetRangeUser(xmin, xmax);
    haxis->GetYaxis()->SetNdivisions(5);
    haxis->GetYaxis()->SetLabelSize(0.12);
    haxis->GetYaxis()->SetLabelOffset(0.01);
    haxis->GetYaxis()->SetTitleSize(0.145);
    haxis->GetYaxis()->SetTitleOffset(0.35);
    haxis->GetXaxis()->SetLabelSize(0.15);
    haxis->GetXaxis()->SetLabelOffset(0.008);
    haxis->GetXaxis()->SetTitleSize(0.2);
    haxis->GetXaxis()->SetTitleOffset(0.8);
    haxis->SetTitle("");
    haxis->GetYaxis()->SetTitle("Ratio");
    haxis->GetXaxis()->SetTitle((unit == "") ? xtitle.Data() : (xtitle + " (" + unit + ")").Data());
    haxis->Draw();

    TLine* line = new TLine(xmin, 1., xmax, 1.);
    line->SetLineWidth(2);
    line->SetLineStyle(kDashed);
    line->SetLineColor(kBlack);
    line->Draw("same");
  }

  pad1->cd();
  Mu2e_lumi(false, plot.npot_, plot.livetime_);

  return c;
}

//------------------------------------------------------------------------------------------------------------
TCanvas* plot_fit_frame(RooPlot* frame, RooRealVar& obs, TString xtitle, TString ytitle, TString data, TString pdf,
                        const double npot = 3.6e20, const double livetime = 3.e7) {
  if(!frame) return nullptr;
  TGaxis::SetExponentOffset(-0.05, 0.01, "Y"); // ensure the exponent doesn't hit label info

  TCanvas* c = new TCanvas("c_model", "c_model", 1200, 1000);
  TPad* pad1 = new TPad("pad1", "pad1", 0., 0.3, 1., 1.0); pad1->Draw();
  TPad* pad2 = new TPad("pad2", "pad2", 0., 0.0, 1., 0.3); pad2->Draw();
  pad1->SetBottomMargin(0.03); pad1->SetTopMargin(0.10); pad1->SetLeftMargin(0.12); pad1->SetRightMargin(0.05);
  pad2->SetBottomMargin(0.35); pad2->SetTopMargin(0.03); pad2->SetLeftMargin(pad1->GetLeftMargin()); pad2->SetRightMargin(pad1->GetRightMargin());
  pad1->SetFillColor(0); pad1->SetTickx(0); pad1->SetTicky(0);
  pad2->SetFillColor(0); pad2->SetTickx(0); pad2->SetTicky(0);

  pad1->cd();
  frame->SetTitle("");
  frame->GetYaxis()->SetTitle((ytitle == "") ? Form("Events / ( %.1f MeV/c)", frame->GetXaxis()->GetBinWidth(1)) : ytitle);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetLabelOffset(0.01);
  frame->GetYaxis()->SetTitleSize(0.06);
  frame->GetYaxis()->SetTitleOffset(0.95);
  frame->GetYaxis()->SetMaxDigits(3);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetLabelOffset(0.007);
  frame->GetXaxis()->SetTitleSize(0.07);
  frame->GetXaxis()->SetTitleOffset(0.80);
  frame->SetTitle("");
  frame->SetXTitle("");
  frame->Draw();

  pad2->cd();
  auto lower_frame = obs.frame();

  auto hpull = frame->pullHist(data.Data(), pdf.Data());
  if(hpull) {
    hpull->SetName("ratio_fnc");
    // hpull->SetLineColor(kBlack);
    // hpull->SetLineWidth(2);
    // hpull->SetMarkerStyle(20);
    // hpull->SetMarkerSize(0.8);
    lower_frame->addPlotable(hpull,"PE1");
  } else {
    cout << __func__ << ": No pull histogram found with data = " << data.Data() << " pdf = " << pdf.Data() << endl;
  }

  lower_frame->GetYaxis()->SetRangeUser(-3.9,3.9);
  lower_frame->GetYaxis()->SetNdivisions(5);
  lower_frame->GetYaxis()->SetLabelSize(0.125);
  lower_frame->GetYaxis()->SetLabelOffset(0.01);
  lower_frame->GetYaxis()->SetTitleSize(0.145);
  lower_frame->GetYaxis()->SetTitleOffset(0.35);
  lower_frame->GetXaxis()->SetLabelSize(0.15);
  lower_frame->GetXaxis()->SetLabelOffset(0.008);
  lower_frame->GetXaxis()->SetTitleSize(0.2);
  lower_frame->GetXaxis()->SetTitleOffset(0.8);
  lower_frame->SetTitle("");
  lower_frame->GetYaxis()->SetTitle("#frac{N_{data} - N_{fit}}{#sigma_{data}}");
  lower_frame->GetXaxis()->SetTitle(xtitle.Data());
  lower_frame->Draw();

  TLine* line = new TLine(lower_frame->GetXaxis()->GetXmin(), 0., lower_frame->GetXaxis()->GetXmax(), 0.);
  line->SetLineWidth(2);
  line->SetLineStyle(kDashed);
  line->SetLineColor(kBlack);
  line->Draw("same");

  pad1->cd();
  auto logo = Mu2e_lumi(false, npot, livetime);
  logo->Draw();

  return c;
}

//------------------------------------------------------------------------------------------------------------
void print_pdf(RooAbsPdf* pdf) {
  printf("Print PDF %s\n", pdf->GetName());
  auto vars = RooArgList(*(pdf->getVariables()));
  for(auto obj : vars) {
    auto var = (RooRealVar*) obj;
    printf(" %40s: %10f +- %10f\n", var->GetName(), var->getVal(), var->getError());
  }

}


//------------------------------------------------------------------------------------------------------------
double closure_integral(double x1, double x2 = 1.) {
  if(x1 >= x2) return 0.;
  if(x1 < 0. || x2 > 1.) return 0.;
  const double y1 = 20./3.*pow(x1,6)-24.*pow(x1,5)+35.*pow(x1,4)-80./3.*pow(x1,3)+10.*pow(x1,2);
  const double y2 = 20./3.*pow(x2,6)-24.*pow(x2,5)+35.*pow(x2,4)-80./3.*pow(x2,3)+10.*pow(x2,2);
  return (y2 - y1);
}

#endif
