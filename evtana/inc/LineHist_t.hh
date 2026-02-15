//
// Line histograms
// Michael MacKenzie (2026)

#ifndef RUN1BEVTANA_LINEHIST_T_HH
#define RUN1BEVTANA_LINEHIST_T_HH

// ROOT includes
#include "TH1.h"
#include "TH2.h"

namespace Run1BEvtAna {
  struct LineHist_t {
    TH1*    fT0;
    TH1*    fT0Err;
    TH1*    fD0;
    TH1*    fChi2NDof;
    TH1*    fFitCons[2];
    TH1*    fTanDip;
    TH1*    fNActive;
    TH1*    fClusterE;
    TH1*    fDt;
    TH1*    fTrackID;
    TH1*    fExlTrackID; //flagging what is exclusively cut by a given ID bit

    //MC truth information
    TH1*    fMCPdg[2];
    TH1*    fMCStrawHits;
    TH1*    fMCGoodHits;
    TH1*    fMCTrajectory;
    TH1*    fMCSimProc;
  };
}
#endif
