#ifndef __CONVANA_ANALYSIS_TYPES__
#define __CONVANA_ANALYSIS_TYPES__

struct pdf_info {
  RooAbsPdf* pdf_ = nullptr;
  TH1* hist_ = nullptr;
  RooRealVar* obs_ = nullptr; //observable variable
  RooRealVar* norm_ = nullptr; //normalization variable
  double rate_ = 0.;
  int color_ = kRed;
  TString name_ = "";
  TString title_ = "";
  int sys_ = -1;
};

struct DatasetInfo_t {
  double ngen_; //for beam processes this is N(gen evt), for cosmics this is gen livetime
  Long64_t ndigi_; //N(events) expected in the digi dataset for missing data corrections
  double emin_;
  double emax_;
  double theory_; //normalization info, without lumi/livetime
  TString name_; //ntuple dataset name
  TString dsname_; //input digi dataset name in SAM
  DatasetInfo_t(double ngen = 1., Long64_t ndigi = 0, double emin = 0., double emax = 1., double theory = 1.,
                TString name = "default", TString dsname = "none") : ngen_(ngen), ndigi_(ndigi), emin_(emin), emax_(emax),
                                                                     theory_(theory), name_(name), dsname_(dsname) {}
  double norm(const double total_rate = 1.) {
    if(ngen_ <= 0.) {
      cout << "DatasetInfo_t::" << __func__ << ": Bad norm input for " << name_.Data() << std::endl;
      return 1.;
    }
    double val(1.);
    if(emin_ < emax_) val *= (emax_ - emin_); //multiply the flat PDF weight
    val /= ngen_; //divide out the generated statistics
    val *= theory_; //apply cross sections/stopping rates/etc.
    val *= total_rate;
    return val;
  }
};

struct plot_t {
  TString hist_      = ""   ;
  TString type_      = ""   ;
  int     selection_ =  0   ;
  double xmin_       =  1.  ;
  double xmax_       = -1.  ;
  double ymin_       =  1.  ;
  double ymax_       = -1.  ;
  int    rebin_      =  1   ;
  bool   logx_       = false;
  bool   logy_       = false;
  TString unit_      = ""   ;
  TString xtitle_    = ""   ;
  TString ytitle_    = ""   ;
  TString title_     = ""   ;
  int     sys_up_    = -1   ;
  int     sys_down_  = -1   ;
  plot_t(TString hist, TString type, int selection,
         int rebin = 1, double xmin = 1., double xmax = -1., double ymin = 1., double ymax = -1.,
         bool logy = false, bool logx = false,
         TString xtitle = "", TString unit = "", TString ytitle = "", TString title = "") :
    hist_(hist), type_(type), selection_(selection), rebin_(rebin), xmin_(xmin), xmax_(xmax),
    ymin_(ymin), ymax_(ymax), logx_(logx), logy_(logy), xtitle_(xtitle), unit_(unit), ytitle_(ytitle), title_(title)
  {}
};


using pdf_pair = std::pair<RooAbsPdf*, double>;
using pdf_pairs = std::vector<pdf_pair>;
using pdf_infos = std::vector<pdf_info>;

#endif
