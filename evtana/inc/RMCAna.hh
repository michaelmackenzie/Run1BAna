//
// RMCAna: RMC analysis ntupling
// Michael MacKenzie (2025)

#ifndef RUN1BEVTANA_RMCANA_HH
#define RUN1BEVTANA_RMCANA_HH

// standard includes

// ROOT includes

// Mu2e Offline includes

// Mu2e EventNtuple includes

// local includes
#include "Run1BAna/evtana/inc/Run1BEvtAna.hh"

using namespace mu2e;
namespace Run1BEvtAna {
  class RMCAna : public Run1BEvtAna {
  public:
    RMCAna(int verbose = 0);
    ~RMCAna() {};

    void InitHistSelections();
    bool ProcessEvent();
    void InitializeEvent();

    int InitializeInput();
    int InitializeOutput();

    TString OutputFileName() { return "rmcana_" + name_ + ".root"; }

  };
}

#endif
