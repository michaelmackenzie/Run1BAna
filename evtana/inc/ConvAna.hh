//
// ConvAna: Conversion search analysis ntupling/histogramming
// Michael MacKenzie (2025)

#ifndef RUN1BEVTANA_CONVANA_HH
#define RUN1BEVTANA_CONVANA_HH

// standard includes

// ROOT includes

// Mu2e Offline includes

// Mu2e EventNtuple includes

// local includes
#include "Run1BAna/evtana/inc/Run1BEvtAna.hh"
#include "Run1BAna/evtana/inc/SysHist_t.hh"
#include "Run1BAna/evtana/inc/Systematics.hh"

using namespace mu2e;
namespace Run1BEvtAna {
  class ConvAna : public Run1BEvtAna {
  public:
    enum {kCRVOffset = 1000, kTimeOffset = 2000};
    ConvAna(int verbose = 0);
    ~ConvAna() {};

    void InitHistSelections();
    void BookSystematicHist(SysHist_t* Hist, const char* Folder);
    void BookHistograms(TDirectory* dir);
    bool ProcessEvent();
    void InitializeEvent();

    int InitializeInput();
    int InitializeOutput();

    TString OutputFileName() { return "convana_" + name_ + ".root"; }

    Bool_t             fill_verbose_sys_ = false        ; // add additional info with each systematic

    SysHist_t*         sys_hists_[kMaxHists]            ; // systematic histograms
    TDirectory*        sys_dirs_ [kMaxHists]            ;
    Systematics        systematics_                     ; // systematic information
  };
}

#endif
