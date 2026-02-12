//
// CRV histograms
// Michael MacKenzie (2025)

#ifndef RUN1BEVTANA_CRVHIST_T_HH
#define RUN1BEVTANA_CRVHIST_T_HH

// ROOT includes
#include "TH1.h"
#include "TH2.h"

namespace Run1BEvtAna {
  struct CRVHist_t {
    TH1*    fSector;
    TH1*    fFirstBar;                           // bar # of the first pulse
    TH1*    fNPulses;
    TH1*    fNPe;                                // N(PE) - apparently, the sum
    TH1*    fNPePP;                              // N(PE) per pulse
    TH1*    fStartTime;
    TH1*    fEndTime;
    TH1*    fWidth;
    TH2*    fXVsZ;
    TH2*    fYVsZ;
    TH1*    fCorrTime;
    TH1*    fCorrTimeProp;
    TH1*    fCorrTimeToF;
    TH1*    fApproxTimeST;
    TH1*    fApproxTimeCalo;
    TH1*    fApproxTimeExtrap;
    TH1*    fApproxTimeSTToFront;
    TH1*    fApproxTimeCaloToFront;
    TH1*    fApproxTimeExtrapToFront;
    TH1*    fBarsOneEnd;
    TH1*    fCrvPropdT;
    TH1*    fNSectors;
    TH1*    fBarsTwoEnd;
    TH1*    fNDiffLSectors;
    TH1*    fStubSlope;
    TH1*    fStubSlopeChi2;
    TH1*    fStubSlopeDelta;
    TH1*    fStubQN;
    TH1*    fStubSlopeMCProduct;
  };
}
#endif
