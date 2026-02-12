//
// Calo cluster information
// Michael MacKenzie (2025)

#ifndef RUN1BEVTANA_CALOLUSTER_T_HH
#define RUN1BEVTANA_CALOLUSTER_T_HH

// ROOT includes
#include "Rtypes.h"

// EventNtuple includes
#include "EventNtuple/inc/CaloClusterInfo.hh"
#include "EventNtuple/inc/CaloHitInfo.hh"

namespace Run1BEvtAna {
  struct CaloCluster_t {
    mu2e::CaloClusterInfo* cluster_;


    //-------------------------------------------------
    // Accessors


    //-------------------------------------------------
    // Additional functions


    void Reset() {
      cluster_ = nullptr;
    }

    CaloCluster_t() { Reset(); }
  };
}
#endif
