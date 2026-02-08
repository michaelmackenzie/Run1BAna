#ifndef __RUN1BANA_ANALYSIS_SIGNALMODEL__
#define __RUN1BANA_ANALYSIS_SIGNALMODEL__

#include "Run1BAna/analysis/defaults.C"
#include "Run1BAna/tools/utilities.C"

//---------------------------------------------------------------------------------------------------------------------------
// Process options: mumem, mumep, fele, fpos
int fit_resolution(TString process, const int selection) {
  // Retrieve the input file
  TString stub;
  if     (process == "mumem") stub = "cele0";
  else if(process == "mumep") stub = "cpos0";
  else if(process == "fele" ) stub = "fele0";
  else if(process == "fpos" ) stub = "fpos0";
  else {
    cout << __func__ << ": Unknown process " << process.Data() << endl;
  }
  TFile* f = TFile::Open(Form("%sConvAna.%s.%sb1s5r0000.m%i.hist", hist_path_, hist_func_, hist_mode_, stub.Data()), "READ");
  if(!f) return -1;

  // Retrieve the input histogram
  TH1* h = (TH1*) f->Get(Form("Ana/ConvAna_ConvAna/Hist/trk_%i/dP", selection));
  if(!h) {
    cout << __func__ << ": Input histogram for selection " << selection << " not found in file " << f->GetName() << endl;
    f->Close();
    return 1;
  }
  h = (TH1*) h->Clone(Form("%s_%i", process.Data(), selection));
  h->SetDirectory(0);
  f->Close();

  //----------------------------------------------
  // Construct the fit objects

  // Create the observable
  const float xmin(-4), xmax(4.);
  RooRealVar obs(Form("obs_%i", selection), "#deltap", (xmin+xmax)/2., xmin, xmax, "MeV/c");
  obs.setBins(h->FindBin(xmax-1.e-4) - h->FindBin(xmin+1.e-4) + 1);

  // Create the histogram data
  h = trim_hist(h, xmin, xmax);
  RooDataHist data_hist("data_hist", "Data hist", obs, h);

  // Create the PDF
  RooRealVar* res_mean     = new RooRealVar("res_mean"  , "mean"  , 0.00, -1., 1.);
  RooRealVar* res_sigma    = new RooRealVar("res_sigma" , "sigma" , 0.15, 0., 5.);
  RooRealVar* res_alpha1   = new RooRealVar("res_alpha1", "alpha1", 1.2, 0.1, 10.);
  RooRealVar* res_alpha2   = new RooRealVar("res_alpha2", "alpha2", 1.9, 0.1, 10.);
  RooRealVar* res_n1       = new RooRealVar("res_n1"    , "enne1" , 2.7, 0.1, 30.);
  RooRealVar* res_n2       = new RooRealVar("res_n2"    , "enne2" , 3.8, 0.1, 30.);
  RooAbsPdf* res_pdf       = new RooCrystalBall("res_pdf", "Resolution PDF", obs, *res_mean, *res_sigma, *res_alpha1, *res_n1, *res_alpha2, *res_n2);

  //----------------------------------------------
  // Perform the fit

  res_pdf->fitTo(data_hist);

  //----------------------------------------------
  // Plot the results

  auto frame = obs.frame();
  data_hist.plotOn(frame, RooFit::Name("data"));
  res_pdf->plotOn(frame, RooFit::Name("pdf"));
  auto c = plot_fit_frame(frame, obs, "p(Reco) - p(True) (MeV/c)", "", "data", "pdf");
  if(!c) return 10;

  const char* figdir = "figures/resolution";
  gSystem->Exec(Form("[ ! -d %s ] && mkdir -p %s", figdir, figdir));

  c->SaveAs(Form("%s/res_fit_%s_%i.png", figdir, process.Data(), selection));
  const float ymax = h->GetMaximum();
  frame->GetYaxis()->SetRangeUser(1.e-5*ymax, 5.*ymax);
  TPad* pad1 = (TPad*) c->GetPrimitive("pad1");
  if(pad1) pad1->SetLogy();
  c->SaveAs(Form("%s/res_fit_%s_%i_log.png", figdir, process.Data(), selection));
  delete frame;
  delete c;

  print_pdf(res_pdf);
  return 0;
}
#endif
