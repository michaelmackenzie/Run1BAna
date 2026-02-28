//
// Cluster histograms
// Michael MacKenzie (2025)

#ifndef RUN1BEVTANA_CLUSTERHIST_T_HH
#define RUN1BEVTANA_CLUSTERHIST_T_HH

// ROOT includes
#include "TH1.h"
#include "TH2.h"

namespace Run1BEvtAna {
  struct ClusterHist_t {
    TH1*    fDiskID;
    TH1*    fEnergy;
    TH1*    fT0;
    TH1*    fRow;
    TH1*    fCol;
    TH1*    fX;
    TH1*    fY;
    TH1*    fZ;
    TH1*    fR;
    TH1*    fNCr0; // all clustered
    TH1*    fNCr1; // above 1MeV
    TH1*    fYMean;
    TH1*    fZMean;
    TH1*    fSigY;
    TH1*    fSigZ;
    TH1*    fSigR;
    TH1*    fFrE1;
    TH1*    fFrE2;
    TH1*    fSigE1;
    TH1*    fSigE2;
    TH1*    fTimeRMS;
    TH1*    fMaxR;
    TH1*    fE9OverE;
    TH1*    fE25OverE;
    TH1*    fRingEOverE;
    TH1*    fRingEOverE1;
    TH1*    fOutRingE;
    TH1*    fOutRingEOverE;

    TH1*    fMCSimEDep;
    TH1*    fMCSimMomIn;
    TH1*    fMCSimPdg;
    TH1*    fMCEDep;
    TH1*    fMCTime;
    TH1*    fMC_dE;
    TH1*    fMC_dt;
  };
}
#endif
