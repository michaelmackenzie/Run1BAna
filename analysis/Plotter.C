#ifndef __RUN1BANA_PLOTTER__
#define __RUN1BANA_PLOTTER__
// Handle plotting distributions

#include "Run1BAna/analysis/dataset_info.C"
#include "Run1BAna/analysis/types.C"
#include "Run1BAna/analysis/defaults.C"
#include "Run1BAna/analysis/physics.C"
#include "Run1BAna/tools/utilities.C"

class Plotter {
public:

  //-------------------------------------------------------------------------------------------------------
  // data structs for plotting

  struct data_t {
    double norm_ = 1.; //normalization scale
    double scale_= 1.; //additional scale factor
    int    set_offset_ = 0; //offset for distribution shapes
    bool   use_set_norm_ = false; // get normalization from nominal set, use shape from offset set
    TFile* f_ = nullptr;
    TString name_ = "default";
    TString label_ = "default";
    int     type_  = 0;
    int     color_ = kRed;
    data_t(double norm = 1., TString name = "default", TString label = "default", int type = 0,
           int color = kRed, int offset = 0, double scale = 1., bool use_set_norm = false) {
      init(norm, name, label, type, color, offset, scale, use_set_norm);
    }
    void init(double norm, TString name, TString label, int type, int color, int offset, double scale, bool use_set_norm = false) {
      norm_ = norm; name_ = name; label_ = label; type_ = type; color_ = color;
      set_offset_ = offset; scale_ = scale; use_set_norm_ = use_set_norm;
    }
  };

  //-------------------------------------------------------------------------------------------------------
  // Enums

  enum {kRatio, kDifference};

  //-------------------------------------------------------------------------------------------------------
  // Fields
  vector<data_t> data_                   ; // input histogram data
  int debug_        = 0                  ; // for debug printout
  TString figdir_   = "figures/plots"    ; // figure directory for plots
  int stack_signal_ = 0                  ; // include the signal model in the model vs. data comparisons
  bool use_offsets_ = true               ; // use control region offsets and scales
  int comp_plot_    = kRatio             ; // ratio or difference plot
  double min_ratio_ = 0.5                ; // ratio plot range
  double max_ratio_ = 1.5                ;
  bool ad_hoc_sys_  = true               ; // ad-hoc systematic band
  TString dir_path_ = "Run1BAna"         ; // directory path in the ROOT file
  vector<TString> bkgs_ = {"dio_90", "dio_80"}; // backgrounds included in the plotting

  //-------------------------------------------------------------------------------------------------------
  Plotter() {}

  //-------------------------------------------------------------------------------------------------------
  Long64_t get_nsampled(const data_t& input) {
    if(!input.f_) return 1;
    TString hist_name = Form("%s/data/norm", dir_path_.Data());
    TH1* h = (TH1*) input.f_->Get(hist_name.Data());
    if(!h) {
      std::cout << "Unable to find " << hist_name.Data() << " in input " << input.name_.Data() << std::endl;
      return 1;
    }
    return h->GetEntries();
  }

  //-------------------------------------------------------------------------------------------------------
  vector<TH1*> get_histograms(TString hist, TString type, int selection, int data_type = -10, TString name_tag = "") {
    map<TString, TH1*> hist_map;
    vector<TString> labels; //force a fixed order
    for(auto& input : data_) {
      if(debug_ > 5) printf("%s: Checking input %s (%s)\n", __func__, input.name_.Data(), input.label_.Data());
      if(data_type >= -1 && input.type_ != data_type) continue;
      if(name_tag != "" && !input.name_.Contains(name_tag)) continue;
      if(!input.f_) {
        cout << __func__ << ": Input file for " << input.name_.Data() << " (" << input.label_.Data() << ") is not open\n";
        continue;
      }
      int set_offset = input.set_offset_;
      if(set_offset > 0 && selection > set_offset && (selection / set_offset) % 2 == 1) set_offset = 0; // in the control region
      const int hist_set = (use_offsets_) ? selection + set_offset : selection;
      if(debug_ > 5) printf("%s: Retrieving input %s (%s) from %s_%i\n", __func__, input.name_.Data(), input.label_.Data(), type.Data(), hist_set);
      TH1* h = (TH1*) input.f_->Get(Form("%s/%s_%i/%s", dir_path_.Data(), type.Data(), hist_set, hist.Data()));
      if(!h) {
        cout << __func__ << ": Histogram " << hist.Data() << "/" << type.Data() << "/" << hist_set
             << " for " << input.name_.Data() << " (" << input.label_.Data() << ") is not found\n";
        continue;
      }
      h = (TH1*) h->Clone("tmp");
      h->SetDirectory(0);
      double scale = input.norm_*(((set_offset > 0 && use_offsets_) || input.type_ < 0) ? input.scale_ : 1.);
      if(data_type != 0) {
        if(input.name_.Contains("cosmic")) scale *= livetime_;
        else                               scale *= npot_    ;
      }
      if(input.use_set_norm_ && set_offset > 0) {
        double offset_norm = 1.;
        TH1* h_nom = (TH1*) input.f_->Get(Form("%s/%s_%i/%s", dir_path_.Data(), type.Data(), selection, hist.Data()));
        if(!h_nom) {
          cout << __func__ << ": Histogram " << hist.Data() << "/" << type.Data() << "/" << selection
               << " for " << input.name_.Data() << " (" << input.label_.Data() << ") nominal set not found\n";
        } else {
          const double n_offset  = h->Integral(0, h->GetNbinsX()+1);
          const double n_nominal = h_nom->Integral(0, h_nom->GetNbinsX()+1);
          offset_norm = (n_offset > 0.) ? n_nominal / n_offset : 0.;
        }
        if(debug_ > 6) printf("%s: Input %s (%s) using nominal / offset scale = %.3g\n", __func__, input.name_.Data(), input.label_.Data(), offset_norm);
        scale *= offset_norm;
      }
      if(debug_ > 6) printf("%s: Scaling input %s (%s) by %.3e\n", __func__, input.name_.Data(), input.label_.Data(), scale);
      h->Scale(scale);
      if(debug_ > 5) printf("%s: Input %s (%s) has integral %.3e\n", __func__, input.name_.Data(), input.label_.Data(), h->Integral());

      if(hist_map.count(input.label_)) {
        hist_map[input.label_]->Add(h);
        delete h;
      } else {
        if(debug_ > 4) printf("%s: Found new label %s\n", __func__, input.label_.Data());
        labels.push_back(input.label_);
        h->SetName(Form("h_%s_%s_%i_%s", hist.Data(), type.Data(), selection, input.name_.Data()));
        h->SetTitle(input.label_);
        h->SetLineColor(input.color_);
        h->SetFillColor(input.color_);
        h->SetLineWidth(2);
        if(input.type_ == 1) {
          h->SetFillStyle(1001);
        } else if(input.type_ == -1) {
          h->SetFillStyle(3005);
        } else {
          h->SetFillStyle(0);
        }
        hist_map[input.label_] = h;
      }
    }
    vector<TH1*> histograms;
    for(auto label : labels) {
      if(debug_ > 4) printf("%s: Adding label %s: Integral = %.3e\n", __func__, label.Data(), hist_map[label]->Integral());
      histograms.push_back(hist_map[label]);
    }
    return histograms;
  }

  //-------------------------------------------------------------------------------------------------------
  TH1* get_histogram(TString hist, TString type, int selection, int data_type, TString name_tag = "") {
    auto hists = get_histograms(hist, type, selection, data_type, name_tag);
    if(hists.empty()) return nullptr;
    if(!hists[0]) return nullptr;
    TH1* h = (TH1*) hists[0]->Clone(Form("h_%s_%s_%i_%i%s", hist.Data(), type.Data(), selection, data_type, (name_tag == "") ? "" : ("_" + name_tag).Data()));
    for(size_t index = 1; index < hists.size(); ++index) {
      if(!hists[index]) continue;
      h->Add(hists[index]);
    }
    return h;
  }

  //-------------------------------------------------------------------------------------------------------
  TCanvas* plot_systematic(plot_t plot) {
    TString hist       = plot.hist_      ;
    TString type       = plot.type_      ;
    int     selection  = plot.selection_ ;
    int     sys_up     = plot.sys_up_    ;
    int     sys_down   = plot.sys_down_  ;
    double  xmin       = plot.xmin_      ;
    double  xmax       = plot.xmax_      ;
    double  ymin       = plot.ymin_      ;
    double  ymax       = plot.ymax_      ;
    int     rebin      = plot.rebin_     ;
    bool    logx       = plot.logx_      ;
    bool    logy       = plot.logy_      ;
    TString unit       = plot.unit_      ;
    TString xtitle     = plot.xtitle_    ;
    TString ytitle     = plot.ytitle_    ;
    TString title      = plot.title_     ;

    if(debug_ > 0) printf("%s: Plotting systematic %s/%s/%i: up = %i, down = %i\n", __func__, hist.Data(), type.Data(), selection,
                          sys_up, sys_down);
    if(debug_ > 1) printf(" Inputs: xmin = %.1f, xmax = %.1f, unit = %s, xtitle = %s\n",
                          xmin, xmax, unit.Data(), xtitle.Data());

    if(data_.empty()) {
      printf("%s: Data not yet initialized, first run \"Plotter::init\"\n", __func__);
      return nullptr;
    }

    if(sys_up < 0) {
      printf("%s: Systematic up index %i is undefined\n", __func__, sys_up);
    }

    // Retrieve the nominal distributions
    vector<TH1*> backgrounds = get_histograms(hist, type, selection,  1);
    vector<TH1*> signals     = get_histograms(hist, type, selection, -1);
    vector<TH1*> datas       = get_histograms(hist, type, selection,  0);
    if(debug_ > 0) cout << __func__ << ": Backgrounds size = " << backgrounds.size() << endl;
    if(debug_ > 0) cout << __func__ << ": Signals size = " << signals.size() << endl;
    if(debug_ > 0) cout << __func__ << ": Datas size = " << datas.size() << endl;

    if(backgrounds.empty()) {
      printf("%s: Background model not found for %s/%s/%i: up = %i, down = %i\n", __func__, hist.Data(), type.Data(), selection,
             sys_up, sys_down);
      return nullptr;
    }

    // Retrieve the systematic variations
    vector<TH1*> bkgs_up   = get_histograms(hist+Form("_%i", sys_up), "sys", selection,  1);
    vector<TH1*> sigs_up   = get_histograms(hist+Form("_%i", sys_up), "sys", selection, -1);
    vector<TH1*> bkgs_down = (sys_down > 0) ? get_histograms(hist+Form("_%i", sys_down), "sys", selection,  1) : vector<TH1*>{};
    vector<TH1*> sigs_down = (sys_down > 0) ? get_histograms(hist+Form("_%i", sys_down), "sys", selection, -1) : vector<TH1*>{};

    if(bkgs_up.empty()) {
      printf("%s: Up Background model not found for %s/%s/%i: up = %i, down = %i\n", __func__, hist.Data(), type.Data(), selection,
             sys_up, sys_down);
      return nullptr;
    }

    // Retrieve the data histograms
    TH1* data = (datas.empty()) ? nullptr : (TH1*) datas[0]->Clone(Form("d_%s_%s_%i", hist.Data(), type.Data(), selection));
    for(size_t idata = 1; idata < datas.size(); ++idata) {
      if(rebin > 1) data->Rebin(rebin);
      data->Add(datas[idata]);
    }
    if(data) {
      data->SetLineWidth(2);
      data->SetMarkerStyle(20);
      data->SetMarkerSize(1.0);
      data->SetLineColor(kBlack);
      data->SetMarkerColor(kBlack);
    }

    // Histograms of the total background model
    TH1* bkg_total = (TH1*) backgrounds[0]->Clone(Form("nom_%s_%s_%i", hist.Data(), type.Data(), selection));
    bkg_total->Reset(); // clear to make adding easier
    for(auto bkg : backgrounds) {
      if(rebin > 1) bkg->Rebin(rebin);
      bkg_total->Add(bkg);
    }
    TH1* bkg_only = (TH1*) bkg_total->Clone(Form("bkg_only_%s_%s_%i", hist.Data(), type.Data(), selection));
    for(auto signal : signals) {
      if(rebin > 1) signal->Rebin(rebin);
      if(stack_signal_) bkg_total->Add(signal);
    }

    // Histograms of the total systematic background model
    TH1* bkg_total_up = (TH1*) bkgs_up[0]->Clone(Form("up_%s_%s_%i", hist.Data(), type.Data(), selection));
    bkg_total_up->Reset(); // clear to make adding easier
    for(auto bkg : bkgs_up) {
      if(rebin > 1) bkg->Rebin(rebin);
      bkg_total_up->Add(bkg);
      if(debug_ > 1) printf("Up: Adding process %s: Integral = %10.3f\n", bkg->GetTitle(), bkg->Integral());
    }
    for(auto signal : sigs_up) {
      if(rebin > 1) signal->Rebin(rebin);
      if(stack_signal_) bkg_total_up->Add(signal);
    }
    TH1* bkg_total_down = nullptr;
    if(!bkgs_down.empty()) {
      bkg_total_down = (TH1*) bkgs_down[0]->Clone(Form("down_%s_%s_%i", hist.Data(), type.Data(), selection));
      bkg_total_down->Reset(); // clear to make adding easier
      for(auto bkg : bkgs_down) {
        if(rebin > 1) bkg->Rebin(rebin);
        bkg_total_down->Add(bkg);
      }
      for(auto signal : sigs_down) {
        if(rebin > 1) signal->Rebin(rebin);
        if(stack_signal_) bkg_total_down->Add(signal);
      }
    }

    // TGraph of the value and uncertainty
    auto graph = errors_from_th1(bkg_total, bkg_total_up, bkg_total_down);

    if(debug_ > 0) {
      printf("Background yields: nominal = %.3g, up = %.3g, down = %.3g\n",
             bkg_total->Integral(), bkg_total_up->Integral(), (bkg_total_down) ? bkg_total_down->Integral() : 0.);
    }
    if(debug_ > 2) {
      printf("%s: Printing background bins:\n", __func__);
      printf("Bin :    Nominal         Up       Down         g0        gup      gdown\n");
      for(int ibin = 1; ibin <= bkg_total->GetNbinsX(); ++ibin) {
        printf(" %3i: %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n",
               ibin, bkg_total->GetBinContent(ibin), bkg_total_up->GetBinContent(ibin),
               (bkg_total_down) ? bkg_total_down->GetBinContent(ibin) : 0.,
               graph->GetY()[ibin-1], graph->GetErrorYhigh(ibin-1), graph->GetErrorYlow(ibin-1)
               );
      }
    }

    // Create the canvas
    TCanvas* c = new TCanvas(Form("c_sys_%s_%s_%i_%i", hist.Data(), type.Data(), selection, sys_up), "Canvas", 1200, 1000);
    TPad* pad1 = new TPad("pad1", "pad1", 0., 0.3, 1., 1.0); pad1->Draw();
    TPad* pad2 = new TPad("pad2", "pad2", 0., 0.0, 1., 0.3); pad2->Draw();
    pad1->SetBottomMargin(0.03); pad1->SetTopMargin(0.10); pad1->SetLeftMargin(0.12); pad1->SetRightMargin(0.09);
    pad2->SetBottomMargin(0.35); pad2->SetTopMargin(0.04); pad2->SetLeftMargin(pad1->GetLeftMargin()); pad2->SetRightMargin(pad1->GetRightMargin());
    pad1->SetFillColor(0); pad1->SetTickx(1); pad1->SetTicky(1);
    pad2->SetFillColor(0); pad2->SetTickx(1); pad2->SetTicky(1);

    pad1->cd();
    TH1* haxis = bkg_total;
    const double abs_xmin = haxis->GetXaxis()->GetXmin();
    const double abs_xmax = haxis->GetXaxis()->GetXmax();
    if(xmin >= xmax) {
      xmin = abs_xmin;
      xmax = abs_xmax;
    }
    xmin = max(xmin, abs_xmin);
    xmax = min(xmax, abs_xmax);
    if(debug_ > 1) cout << "Axis ranges: xmin = " << xmin << " xmax = " << xmax << endl;

    double max_val = max_in_range(bkg_total, xmin, xmax);
    max_val = max(max_val, max_in_range(bkg_total_up, xmin, xmax));
    if(bkg_total_down) max_val = max(max_val, max_in_range(bkg_total_down, xmin, xmax));
    if(debug_ > 1) cout << "Max val at bkg: " << max_val << endl;
    haxis->SetTitle("");
    haxis->SetXTitle("");
    if(ytitle == "") haxis->SetYTitle(Form("Entries / %.2g %s", haxis->GetBinWidth(1), unit.Data()));
    haxis->GetYaxis()->SetLabelSize(0.05);
    haxis->GetXaxis()->SetLabelSize(0.);
    haxis->GetYaxis()->SetTitleSize(0.06);
    haxis->GetYaxis()->SetTitleOffset(0.9);
    haxis->Draw();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // Draw the model
    bkg_total->SetLineWidth(2);
    bkg_total->SetLineColor(kRed);
    bkg_total->SetFillColor(0);
    bkg_total->Draw("hist");

    if(stack_signal_) {
      bkg_only->SetLineWidth(2);
      bkg_only->SetLineColor(kRed);
      bkg_only->SetFillColor(0);
      bkg_only->SetLineStyle(kDashed);
      bkg_only->Draw("hist same");
    }

    // Draw the uncertainty
    graph->SetFillColor(kRed);
    graph->SetFillStyle(3001);
    graph->Draw("E2");

    // Draw the signals
    for(auto signal : signals) {
      signal->Draw("hist same");
      max_val = max(max_val, max_in_range(signal, xmin, xmax));
      if(debug_ > 1) cout << "Max val at signal: " << max_val << endl;
    }

    // Draw the data
    if(data) {
      data->Draw("EX0 same");
      max_val = max(max_val, max_in_range(data, xmin, xmax, true));
      if(debug_ > 1) cout << "Max val at data: " << max_val << endl;
    }

    if(ymin > ymax) {
      if(logy) { ymin = 1.e-5*max_val; ymax = 50.*max_val;}
      else     { ymin = 1.e-5*max_val; ymax = 1.3*max_val;}
    }
    if(debug_ > 1) cout << "Axis ranges: ymin = " << ymin << " ymax = " << ymax << endl;
    haxis->GetYaxis()->SetRangeUser(ymin, ymax);
    if(xmin < xmax) haxis->GetXaxis()->SetRangeUser(xmin+1.e-6, xmax-1.e-6);
    if(logy) pad1->SetLogy();
    if(logx) pad1->SetLogx();

    // Add a legend
    TLegend* leg = new TLegend(pad1->GetLeftMargin()+0.03, 1. - pad1->GetTopMargin()-0.25, 1.-pad1->GetRightMargin()-0.02, 1. - pad1->GetTopMargin()-0.06);
    leg->SetNColumns(3);
    leg->SetLineWidth(0); leg->SetLineColor(0); leg->SetFillColor(0); leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    if(data) leg->AddEntry(data, "Data", "PL");
    for(auto signal : signals) {
      leg->AddEntry(signal, signal->GetTitle(), "F");
    }
    if(stack_signal_) {
      leg->AddEntry(bkg_only , "Background", "L");
      leg->AddEntry(bkg_total, "Total", "L");
    } else {
      leg->AddEntry(bkg_total, "Background", "L");
    }
    leg->AddEntry(graph, "Systematic", "F");

    leg->Draw();

    // Add a ratio plot
    pad2->cd();
    TH1* haxis_r = (TH1*) haxis->Clone(Form("axis_r_%s_%s_%i", hist.Data(), type.Data(), selection));
    haxis_r->Reset();
    haxis_r->SetFillColor(0); haxis_r->SetLineWidth(0); haxis_r->SetLineColor(0);
    haxis_r->Draw("hist");
    auto graph_r = (TGraphAsymmErrors*) graph->Clone(Form("g_r_%s_%s_%i_%i", hist.Data(), type.Data(), selection, sys_up));
    double max_r(0.), min_r(1.e10);
    for(int ipoint = 0; ipoint < graph_r->GetN(); ++ipoint) {
      const double val = bkg_total->GetBinContent(ipoint+1);
      if(val <= 0.) {
        graph_r->SetPointY     (ipoint, -10.);
        graph_r->SetPointEYhigh(ipoint, 0.);
        graph_r->SetPointEYlow (ipoint, 0.);
        continue;
      }
      graph_r->SetPointY(ipoint, 1.);
      const double up = bkg_total_up->GetBinContent(ipoint+1) / val;
      const double down = (bkg_total_down) ? bkg_total_down->GetBinContent(ipoint+1) / val : 1.;
      graph_r->SetPointEYhigh(ipoint, max(0., max(up-1., down-1.)));
      graph_r->SetPointEYlow (ipoint, max(0., max(1.-up, 1.-down)));
      max_r = max(max_r, graph_r->GetErrorYhigh(ipoint) + 1.);
      min_r = min(min_r, 1. - graph_r->GetErrorYlow (ipoint));
    }
    graph_r->Draw("E2");
    haxis_r->GetYaxis()->SetRangeUser(min_r - 0.05*(max_r - min_r), max_r + 0.05*(max_r - min_r));
    TLine* line = new TLine(xmin, 1., xmax, 1.);
    line->SetLineWidth(2);
    line->SetLineColor(kBlack);
    line->SetLineStyle(kDashed);
    line->Draw("same");

    haxis_r->SetYTitle("Shift / Nominal");
    haxis_r->SetXTitle((unit == "") ? xtitle : xtitle + " (" + unit + ")");
    haxis_r->GetYaxis()->SetNdivisions(507);
    haxis_r->GetYaxis()->SetLabelSize(0.125);
    haxis_r->GetYaxis()->SetLabelOffset(0.01);
    haxis_r->GetYaxis()->SetTitleSize(0.145);
    haxis_r->GetYaxis()->SetTitleOffset(0.35);
    haxis_r->GetXaxis()->SetLabelSize(0.15);
    haxis_r->GetXaxis()->SetLabelOffset(0.008);
    haxis_r->GetXaxis()->SetTitleSize(0.18);
    haxis_r->GetXaxis()->SetTitleOffset(0.8);
    if(logx) pad2->SetLogx();

    // Add logo info
    pad1->cd();
    auto logo = Mu2e_lumi(false, npot_, livetime_);
    logo->Draw();

    // Redraw to get axes on top
    pad1->RedrawAxis();
    pad2->RedrawAxis();

    return c;
  }

  //-------------------------------------------------------------------------------------------------------
  TCanvas* plot_component(plot_t plot, TString process) {
    TString hist       = plot.hist_      ;
    TString type       = plot.type_      ;
    int     selection  = plot.selection_ ;
    double  xmin       = plot.xmin_      ;
    double  xmax       = plot.xmax_      ;
    double  ymin       = plot.ymin_      ;
    double  ymax       = plot.ymax_      ;
    int     rebin      = plot.rebin_     ;
    bool    logx       = plot.logx_      ;
    bool    logy       = plot.logy_      ;
    TString unit       = plot.unit_      ;
    TString xtitle     = plot.xtitle_    ;
    TString ytitle     = plot.ytitle_    ;
    TString title      = plot.title_     ;

    if(debug_ > 0) printf("%s: Plotting distribution %s/%s/%i for component %s\n", __func__, hist.Data(), type.Data(), selection, process.Data());
    if(debug_ > 1) printf(" Inputs: xmin = %.1f, xmax = %.1f, unit = %s, xtitle = %s\n",
                          xmin, xmax, unit.Data(), xtitle.Data());

    if(data_.empty()) {
      printf("%s: Data not yet initialized, first run \"init()\"\n", __func__);
      return nullptr;
    }

    vector<TH1*> hists = get_histograms(hist, type, selection,  -10, process);
    if(debug_ > 0) cout << "Histograms size = " << hists.size() << endl;

    if(hists.empty()) return nullptr;

    // Retrieve the data histograms
    TH1* h = (TH1*) hists[0]->Clone(Form("h_%s_%s_%s_%i", process.Data(), hist.Data(), type.Data(), selection));
    for(size_t ihist = 1; ihist < hists.size(); ++ihist) h->Add(hists[ihist]);
    if(rebin > 1) h->Rebin(rebin);

    // Create the canvas
    TCanvas* c = new TCanvas(Form("c_%s_%s_%s_%i", process.Data(), hist.Data(), type.Data(), selection), "Canvas", 1200, 700);
    c->SetBottomMargin(0.13); c->SetTopMargin(0.10); c->SetLeftMargin(0.12); c->SetRightMargin(0.09);
    c->SetFillColor(0); c->SetTickx(1); c->SetTicky(1);

    const double abs_xmin = h->GetXaxis()->GetXmin();
    const double abs_xmax = h->GetXaxis()->GetXmax();
    if(xmin >= xmax) {
      xmin = abs_xmin;
      xmax = abs_xmax;
    }
    xmin = max(xmin, abs_xmin);
    xmax = min(xmax, abs_xmax);

    const double max_val = max_in_range(h, xmin, xmax);
    if(debug_ > 1) cout << "Max val: " << max_val << endl;
    h->SetTitle("");
    h->SetXTitle("");
    if(ytitle == "") h->SetYTitle(Form("Entries / %.2g %s", h->GetBinWidth(1), unit.Data()));
    else             h->SetYTitle("");
    h->SetXTitle((unit == "") ? xtitle : xtitle + " (" + unit + ")");

    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetLabelOffset(0.005);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(0.9);

    h->Draw("hist");
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    if(ymin > ymax) {
      if(logy) { ymin = 1.e-5*max_val; ymax = 50.*max_val;}
      else     { ymin = 1.e-5*max_val; ymax = 1.3*max_val;}
    }
    if(debug_ > 1) cout << "Axis ranges: ymin = " << ymin << " ymax = " << ymax << endl;
    h->GetYaxis()->SetRangeUser(ymin, ymax);
    if(xmin < xmax) h->GetXaxis()->SetRangeUser(xmin+1.e-6, xmax-1.e-6);
    if(logy) c->SetLogy();
    if(logx) c->SetLogx();

    // Add the logo
    auto logo = Mu2e_lumi(false, npot_, livetime_);
    logo->Draw();

    // Redraw to get axes on top
    c->RedrawAxis();
    c->RedrawAxis();

    return c;
  }


  //-------------------------------------------------------------------------------------------------------
  TCanvas* plot_stack(plot_t plot) {
    TString hist       = plot.hist_      ;
    TString type       = plot.type_      ;
    int     selection  = plot.selection_ ;
    double  xmin       = plot.xmin_      ;
    double  xmax       = plot.xmax_      ;
    double  ymin       = plot.ymin_      ;
    double  ymax       = plot.ymax_      ;
    int     rebin      = plot.rebin_     ;
    bool    logx       = plot.logx_      ;
    bool    logy       = plot.logy_      ;
    TString unit       = plot.unit_      ;
    TString xtitle     = plot.xtitle_    ;
    TString ytitle     = plot.ytitle_    ;
    TString title      = plot.title_     ;

    if(debug_ > 0) printf("%s: Plotting stack %s/%s/%i\n", __func__, hist.Data(), type.Data(), selection);
    if(debug_ > 1) printf(" Inputs: xmin = %.1f, xmax = %.1f, unit = %s, xtitle = %s\n",
                          xmin, xmax, unit.Data(), xtitle.Data());

    if(data_.empty()) {
      printf("%s: Data not yet initialized, first run \"init()\"\n", __func__);
      return nullptr;
    }

    vector<TH1*> backgrounds = get_histograms(hist, type, selection,  1);
    vector<TH1*> signals     = get_histograms(hist, type, selection, -1);
    vector<TH1*> datas       = get_histograms(hist, type, selection,  0);
    if(debug_ > 0) cout << "Backgrounds size = " << backgrounds.size() << endl;
    if(debug_ > 0) cout << "Signals size = " << signals.size() << endl;
    if(debug_ > 0) cout << "Datas size = " << datas.size() << endl;

    if(backgrounds.empty()) return nullptr;

    // Retrieve the data histograms
    TH1* data = (datas.empty()) ? nullptr : (TH1*) datas[0]->Clone(Form("d_%s_%s_%i", hist.Data(), type.Data(), selection));
    for(size_t idata = 1; idata < datas.size(); ++idata) data->Add(datas[idata]);
    if(data) {
      data->SetLineWidth(2);
      data->SetMarkerStyle(20);
      data->SetMarkerSize(1.0);
      data->SetLineColor(kBlack);
      data->SetMarkerColor(kBlack);
    }

    // Create a stack of the processes
    double min_val = 1.e10;
    THStack* stack = new THStack(Form("s_%s_%s_%i", hist.Data(), type.Data(), selection), "Background stack");
    for(auto bkg : backgrounds) {
      if(rebin > 1) bkg->Rebin(rebin);
      stack->Add(bkg);
      min_val = min(min_val, min_in_range(bkg, xmin, xmax, false, 1.e-10));
    }
    for(auto signal : signals) {
      if(rebin > 1) signal->Rebin(rebin);
      if(stack_signal_) stack->Add(signal);
      min_val = min(min_val, min_in_range(signal, xmin, xmax, false, 1.e-10));
    }

    // Create the canvas
    TCanvas* c = new TCanvas(Form("c_%s_%s_%i", hist.Data(), type.Data(), selection), "Canvas", 1200, 1000);
    TPad* pad1 = new TPad("pad1", "pad1", 0., 0.3, 1., 1.0); pad1->Draw();
    TPad* pad2 = new TPad("pad2", "pad2", 0., 0.0, 1., 0.3); pad2->Draw();
    pad1->SetBottomMargin(0.03); pad1->SetTopMargin(0.10); pad1->SetLeftMargin(0.12); pad1->SetRightMargin(0.09);
    pad2->SetBottomMargin(0.35); pad2->SetTopMargin(0.04); pad2->SetLeftMargin(pad1->GetLeftMargin()); pad2->SetRightMargin(pad1->GetRightMargin());
    pad1->SetFillColor(0); pad1->SetTickx(1); pad1->SetTicky(1);
    pad2->SetFillColor(0); pad2->SetTickx(1); pad2->SetTicky(1);

    pad1->cd();
    const char* axis_name = Form("axis_%s_%s_%i", hist.Data(), type.Data(), selection);
    TH1* haxis = (TH1*) ((data) ? data->Clone(axis_name) : stack->GetStack()->Last()->Clone(axis_name));
    haxis->Reset();
    haxis->SetLineWidth(0);
    haxis->SetTitle("");
    haxis->SetXTitle("");
    if(ytitle == "") haxis->SetYTitle(Form("Entries / %.2g %s", haxis->GetBinWidth(1), unit.Data()));
    haxis->GetYaxis()->SetLabelSize(0.05);
    haxis->GetXaxis()->SetLabelSize(0.);
    haxis->GetYaxis()->SetTitleSize(0.06);
    haxis->GetYaxis()->SetTitleOffset(0.9);
    haxis->Draw();
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // get the axis range
    const double abs_xmin = haxis->GetXaxis()->GetXmin();
    const double abs_xmax = haxis->GetXaxis()->GetXmax();
    if(xmin >= xmax) {
      xmin = abs_xmin;
      xmax = abs_xmax;
    } else if(xmax < abs_xmin || xmin > abs_xmax) { // poorly defined bounds
      xmin = abs_xmin;
      xmax = abs_xmax;
    }
    xmin = max(xmin, abs_xmin);
    xmax = min(xmax, abs_xmax);
    if(debug_ > 1) cout << "Axis ranges: xmin = " << xmin << " xmax = " << xmax << endl;
    haxis->GetXaxis()->SetRangeUser(xmin, xmax);

    double max_val = max_in_range((TH1*) stack->GetStack()->Last(), xmin, xmax);
    if(debug_ > 1) cout << "Max val at bkg: " << max_val << endl;

    // Draw the stack
    stack->Draw("hist noclear same");

    // Draw the signals
    for(auto signal : signals) {
      signal->Draw("hist same");
      max_val = max(max_val, max_in_range(signal, xmin, xmax));
      if(debug_ > 1) cout << "Max val at signal: " << max_val << endl;
    }

    // Draw the model statistical uncertainty
    TH1* hbkg_total = (TH1*) stack->GetStack()->Last()->Clone(Form("bkg_%s_%s_%i", hist.Data(), type.Data(), selection));
    hbkg_total->SetFillStyle(3001);
    hbkg_total->SetFillColor(kGray+1);
    hbkg_total->SetLineColor(0);
    hbkg_total->SetLineWidth(0);
    hbkg_total->Draw("same E2");


    // Draw the data
    if(data) {
      if(rebin > 1) data->Rebin(rebin);
      data->Draw("EX0 same");
      max_val = max(max_val, max_in_range(data, xmin, xmax, true));
      min_val = min(min_val, min_in_range(data, xmin, xmax, false, 1.e-10));
      if(debug_ > 1) cout << "Max val at data: " << max_val << " Min val " << min_val << endl;
    }

    if(ymin > ymax) {
      if(logy) {
        if(max_val <= 0.) { ymin = 0.1; ymax = 1.;}
        else {
          min_val = max(min_val, 1.e-5*max_val); // only span at most 5 orders
          const double orders = std::log10(max_val/min_val); // orders spanned
          ymin = std::pow(10, -0.05*orders)*min_val; // ~5% buffer on the bottom
          ymax = std::pow(10,  0.40*orders)*max_val; // ~30% buffer on the bottom
        }
      } else { ymin = 1.e-5*max_val; ymax = 1.3*max_val;}
    }
    if(debug_ > 1) cout << "Axis ranges: ymin = " << ymin << " ymax = " << ymax << endl;
    haxis->GetYaxis()->SetRangeUser(ymin, ymax);
    haxis->GetXaxis()->SetRangeUser(xmin, xmax);
    if(logy) pad1->SetLogy();
    if(logx) pad1->SetLogx();

    // Add a legend
    TLegend* leg = new TLegend(pad1->GetLeftMargin()+0.03, 1. - pad1->GetTopMargin()-0.21, 1.-pad1->GetRightMargin()-0.02, 1. - pad1->GetTopMargin()-0.03);
    leg->SetNColumns(3);
    leg->SetLineWidth(0); leg->SetLineColor(0); leg->SetFillColor(0); leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    if(data) leg->AddEntry(data, "Data", "PL");
    for(auto signal : signals) {
      leg->AddEntry(signal, signal->GetTitle(), "F");
    }
    // reverse the order of backgrounds in the legend
    if(!backgrounds.empty()) {
      for(int ibkg = backgrounds.size() - 1; ibkg >= 0; --ibkg) {
        auto bkg = backgrounds[ibkg];
        leg->AddEntry(bkg, bkg->GetTitle(), "F");
      }
    }

    leg->Draw();

    // Add a ratio plot
    pad2->cd();
    TH1* haxis_r = (TH1*) haxis->Clone(Form("axis_r_%s_%s_%i", hist.Data(), type.Data(), selection));

    TH1* hnum = (data) ? (TH1*) data->Clone(Form("data_r_%s_%s_%i", hist.Data(), type.Data(), selection)) :
      (!signals.empty()) ? (TH1*) signals[0]->Clone(Form("signal_r_%s_%s_%i", hist.Data(), type.Data(), selection)) :
      nullptr;
    if(hnum) {
      haxis_r->Draw();
      TLine* line(nullptr);
      if(comp_plot_ == kDifference) {
        hnum->Add(hbkg_total, -1.);
        const double max_diff = max_in_range(hnum, xmin, xmax, true);
        const double min_diff = min_in_range(hnum, xmin, xmax, true, -1.e10);
        const double diff_buffer = 0.05*(max_diff - min_diff);
        haxis_r->GetYaxis()->SetRangeUser(min_diff - diff_buffer, max_diff + diff_buffer);
        haxis_r->SetYTitle((data) ? (stack_signal_) ? "Data - Total" : "Data - Bkg" : "Signal - Bkg");
      } else {  // default to ratio plot
        hnum->Divide(hbkg_total);
        if(data) {
          haxis_r->GetYaxis()->SetRangeUser(min_ratio_, max_ratio_);
          line = new TLine(xmin, 1., xmax, 1.);
          line->SetLineWidth(2);
          line->SetLineColor(kBlack);
          line->SetLineStyle(kDashed);
          line->Draw("same");
          if(stack_signal_) haxis_r->SetYTitle("Data / Total");
          else              haxis_r->SetYTitle("Data / Bkg");
        } else {
          const double max_ratio = hnum->GetMaximum();
          haxis_r->GetYaxis()->SetRangeUser(0., 1.1*max_ratio);
          haxis_r->SetYTitle("Signal / Bkg");
        }
      }
      hnum->Draw("EX0 same");
      haxis_r->SetXTitle((unit == "") ? xtitle : xtitle + " (" + unit + ")");
      haxis_r->GetYaxis()->SetNdivisions(507);
      haxis_r->GetYaxis()->SetLabelSize(0.125);
      haxis_r->GetYaxis()->SetLabelOffset(0.01);
      haxis_r->GetYaxis()->SetTitleSize(0.145);
      haxis_r->GetYaxis()->SetTitleOffset(0.35);
      haxis_r->GetXaxis()->SetLabelSize(0.15);
      haxis_r->GetXaxis()->SetLabelOffset(0.008);
      haxis_r->GetXaxis()->SetTitleSize(0.18);
      haxis_r->GetXaxis()->SetTitleOffset(0.8);
      if(logx) pad2->SetLogx();

      // Add uncertainty bands
      if(data) {
        const int bin_low = (xmin < xmax) ? hbkg_total->FindBin(xmin + 1.e-6) : 1;
        const int bin_high = (xmin < xmax) ? hbkg_total->FindBin(xmax - 1.e-6) : hbkg_total->GetNbinsX();
        const int sys_bins = bin_high - bin_low + 1;
        double sys_x[sys_bins], sys_y[sys_bins], sys_xerr[sys_bins], sys_yerr[sys_bins], stat_yerr[sys_bins];
        for(int ibin = bin_low; ibin <= bin_high; ++ibin) {
          double val = hbkg_total->GetBinContent(ibin);
          double err = hbkg_total->GetBinError  (ibin);
          const int index = ibin - bin_low;
          sys_x    [index] = hbkg_total->GetBinCenter(ibin);
          sys_xerr [index] = hbkg_total->GetBinWidth(ibin)/2.;
          sys_y    [index] = (comp_plot_ == kDifference) ? 0. : 1.;
          sys_yerr [index] = 0.; // default value
          stat_yerr[index] = 0.;
          if(val <= 0.) continue;
          if(ad_hoc_sys_) {
            double val_dio(0.), val_csm(0.), val_rpc(0.), val_beam(0.);
            for(int iproc = 0; iproc < stack->GetNhists(); ++iproc) {
              TH1* hproc = (TH1*) stack->GetHists()->At(iproc);
              const double proc_val = hproc->GetBinContent(ibin);
              TString proc_name(hproc->GetName());
              if     (proc_name.Contains("dio"   )) val_dio  += proc_val;
              else if(proc_name.Contains("cosmic")) val_csm  += proc_val;
              else if(proc_name.Contains("rpc"))    val_rpc  += proc_val;
              else                                  val_beam += proc_val;
            }
            const double sys = sqrt(pow(0.1*(val_dio+val_rpc+val_beam)/val, 2) + pow(0.2*val_csm/val, 2) + pow(0.025*val_dio/val, 2)
                                    + pow(0.093*val_rpc/val, 2) + pow(0.27*val_rpc/val,2));
            if(comp_plot_ == kDifference) sys_yerr[index] = sys*val;
            else                          sys_yerr[index] = sys;
          }
          if(comp_plot_ != kDifference && val > 0.) err /= val;
          stat_yerr[index] = err;
          sys_yerr [index] = sqrt(pow(err, 2) + pow(sys_yerr[index], 2)); //stat + sys error
        }
        TGraphErrors* gsys = new TGraphErrors(sys_bins, sys_x, sys_y, sys_xerr, sys_yerr);
        gsys->SetFillStyle(3004);
        gsys->SetFillColor(kGray+1);
        gsys->SetLineWidth(0);
        gsys->Draw("E2");
        TGraphErrors* gstat = new TGraphErrors(sys_bins, sys_x, sys_y, sys_xerr, stat_yerr);
        gstat->SetFillStyle(3001);
        gstat->SetFillColor(kGray+1);
        gstat->SetLineWidth(0);
        gstat->Draw("E2");
        if(line) line->Draw("same"); // redraw to put in the front
        hnum->Draw("EX0 same");
      }
    } else {
      // Do just a background plot
      pad1->SetBBoxY1(0.);
      pad1->Draw();
    }

    // Add logo info
    pad1->cd();
    auto logo = Mu2e_lumi(false, npot_, livetime_);
    logo->Draw();

    // Redraw to get axes on top
    pad1->RedrawAxis();
    pad2->RedrawAxis();

    return c;
  }

  //-------------------------------------------------------------------------------------------------------
  TCanvas* plot_stack(TString hist, TString type, int selection,
                      int rebin = 1, double xmin = 1., double xmax = -1., double ymin = 1., double ymax = -1.,
                      bool logy = false, bool logx = false
                      ) {
    plot_t plot(hist, type, selection, rebin, xmin, xmax, ymin, ymax, logy, logx);
    return plot_stack(plot);
  }

  //-------------------------------------------------------------------------------------------------------
  TCanvas* print_component(plot_t plot, TString process) {
    TCanvas* c = plot_component(plot, process);
    if(!c) return nullptr;

    c->SaveAs(Form("%s/comp_%s_%s_%s_%i%s.png", figdir_.Data(), process.Data(), plot.hist_.Data(), plot.type_.Data(), plot.selection_,
                   (plot.logy_) ? "_log" : ""
                   ));
    return c;
  }

  //-------------------------------------------------------------------------------------------------------
  TCanvas* print_stack(plot_t plot) {
    TCanvas* c = plot_stack(plot);
    if(!c) return nullptr;

    c->SaveAs(Form("%s/stack_%s_%s_%i%s.png", figdir_.Data(), plot.hist_.Data(), plot.type_.Data(), plot.selection_,
                   (plot.logy_) ? "_log" : ""
                   ));
    return c;
  }

  //-------------------------------------------------------------------------------------------------------
  TCanvas* print_stack(TString hist, TString type, int selection,
                       int rebin = 1, double xmin = 1., double xmax = -1., double ymin = 1., double ymax = -1.,
                       bool logy = false, bool logx = false
                       ) {
    TCanvas* c = plot_stack(hist, type, selection, rebin, xmin, xmax, ymin, ymax, logy, logx);
    if(!c) return nullptr;

    c->SaveAs(Form("%s/stack_%s_%s_%i%s.png", figdir_.Data(), hist.Data(), type.Data(), selection,
                   (logy) ? "_log" : ""
                   ));
    return c;
  }

  //-------------------------------------------------------------------------------------------------------
  TCanvas* print_systematic(plot_t plot) {
    TCanvas* c = plot_systematic(plot);
    if(!c) return nullptr;

    c->SaveAs(Form("%s/sys_%s_%s_%i_sys_%i%s.png", figdir_.Data(), plot.hist_.Data(), plot.type_.Data(), plot.selection_,
                   plot.sys_up_, (plot.logy_) ? "_log" : ""
                   ));
    return c;
  }

  //-------------------------------------------------------------------------------------------------------
  TCanvas* print_systematic(plot_t plot, int sys_up, int sys_down = -1) {
    plot.sys_up_ = sys_up; plot.sys_down_ = sys_down;
    return print_systematic(plot);
  }

  //-------------------------------------------------------------------------------------------------------
  // Update the signal branching fraction
  void update_signal_br(double signal_br) {
    signal_br_ = signal_br;
    for(auto& input : data_) {
      if(input.type_ != -1) continue;
      TString title = input.label_;
      title = (title.Last('(') > 0) ? title(0, title.Last('(')) : title;
      const double log_br = log10(signal_br_) + 1.e-5; //push slightly up
      const int power = log_br - 1; //take an extra off
      title += Form(" (%.1f x 10^{%i})", signal_br_/pow(10,power), power);
      input.label_ = title;
      input.scale_ = signal_br_;
    }
  }

  //-------------------------------------------------------------------------------------------------------
  // Initialize the inputs
  int init(bool mumem = true, TString dataset = "", TString tag = "") {

    // Check if data is already initialized
    if(!data_.empty()) {
      for(auto& data : data_) if(data.f_) data.f_->Close();
      data_.clear();
    }

    // Check if the figure directory exists
    gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir_.Data(), figdir_.Data()));
    TGaxis::SetMaxDigits(3);
    TGaxis::SetExponentOffset(-0.06, 0.008, "Y");

    // Initialize the input data
    // FIXME: This dataset selection should be given, not hard-coded
    get_dataset_info((bkgs_.empty()) ? "dio_80" : bkgs_[0]); //initialize the dataset information if needed
    TString title; int color;
    set_style((mumem) ? "mumem" : "mumep", title, color);
    if(mumem) data_.push_back(data_t(datasets_["mumem"].norm(), "mumem", title, -1, color, 0, signal_br_));
    else      data_.push_back(data_t(datasets_["mumep"].norm(), "mumep", title, -1, color, 0, signal_br_));
    update_signal_br(signal_br_);
    for(auto name : bkgs_) {
      auto info = datasets_[name];
      set_style(name, title, color);
      int offset(0); double scale(1.); bool use_set_norm(false);
      // if(name == "cosmic") {offset = 1000; scale = 1.; use_set_norm = true;} // just normalize to the nominal yield for now
      data_.push_back(data_t(info.norm(), name, title, 1, color, offset, scale, use_set_norm));
    }

    // Add the data if defined
    if(datasets_.count("data_"+dataset)) {
      cout << "Using Data from dataset " << dataset.Data() << endl;
      data_.push_back(data_t(1., "data_"+dataset, "Data (Sim)", 0, kBlack, 0, 1.));
    }
    init_physics(dataset + "_" + tag);

    // Open each input file
    for(auto& input : data_) {
      TString hist_file = datasets_[input.name_].dsname_;
      hist_file.ReplaceAll("mcs.", "nts.");
      hist_file.ReplaceAll("dig.", "nts.");
      hist_file.ReplaceAll(".art", ".root");
      TString file_name = hist_path_ + hist_file;
      if(debug_ > 0) cout << "Opening file " << file_name.Data() << endl;
      input.f_ = TFile::Open(file_name.Data(), "READ");
      if(!input.f_) {
        cout << "File not found: " << file_name.Data() << "(input name = " << input.name_.Data() << ")" << endl;
        return 1;
      }
      const double nsampled = get_nsampled(input);
      const Long64_t nexpect = datasets_[input.name_].ndigi_;
      if(nexpect != nsampled) {
        const double ratio = nexpect * 1. / nsampled;
        cout << __func__ << ": See " << nsampled << " events but expect " << nexpect << " for file " << file_name.Data()
             << " --> scaling by " << ratio << endl;
        input.norm_ *= ratio;
      }
    }
    return 0;
  }
};

#endif
