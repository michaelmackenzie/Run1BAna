//
// Run1BEvtAna systematic histograms
// Michael MacKenzie (2025)

#ifndef RUN1BEVTANA_SYSHIST_T_HH
#define RUN1BEVTANA_SYSHIST_T_HH

// ROOT includes
#include "TH1.h"
#include "TH2.h"

// local includes
#include "Run1BAna/evtana/inc/GlobalConstants.h"

namespace Run1BEvtAna {

  struct SysHist_t {
    TH1* fObs        [kMaxSystematics];
    TH1* fDeltaObs   [kMaxSystematics];
    TH1* fWeight     [kMaxSystematics];
    TH1* fDeltaWeight[kMaxSystematics];

    SysHist_t() {
      for(int i = 0; i < kMaxSystematics; ++i) {
        fObs        [i] = nullptr;
        fDeltaObs   [i] = nullptr;
        fWeight     [i] = nullptr;
        fDeltaWeight[i] = nullptr;
      }
    }
  };
}
#endif
