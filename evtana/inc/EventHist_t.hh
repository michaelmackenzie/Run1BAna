//
// Event histograms
// Michael MacKenzie (2025)

#ifndef RUN1BEVTANA_EVENTHIST_T_HH
#define RUN1BEVTANA_EVENTHIST_T_HH

// ROOT includes
#include "TH1.h"
#include "TH2.h"

namespace Run1BEvtAna {
  struct EventHist_t {
    TH1*    fInstLumi;
    TH1*    fInstLumiApr; // only filled if apr triggered event
    TH1*    fInstLumiCpr; // only filled if cpr triggered event
    TH1*    fInstLumiAprCpr; // filled if either apr or cpr triggered event
    TH1*    fEventWeight[2];
    TH1*    fNAprHelices;
    TH1*    fNCprHelices;
    TH1*    fNHelices;
    TH1*    fNMatchedHelices;
    TH1*    fNAprTracks;
    TH1*    fNCprTracks;
    TH1*    fNTracks;
    TH1*    fNUeTracks;
    TH1*    fNDmuTracks;
    TH1*    fNUmuTracks;
    TH1*    fNCRVClusters;
    TH1*    fNGoodCRVClusters;
    TH1*    fNonCRVVetoID;
    TH1*    fNGoodTrks;
    TH1*    fNIDTrks;
    TH1*    fNDigis;
    TH1*    fNClusters;

    // Primary process info
    TH1*    fPrimaryCode;
    TH1*    fPrimaryGenE;
  };
}
#endif
