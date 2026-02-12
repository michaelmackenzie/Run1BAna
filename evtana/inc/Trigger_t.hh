//
// Trigger information
// Michael MacKenzie (2025)

#ifndef RUN1BEVTANA_TRIGGER_T_HH
#define RUN1BEVTANA_TRIGGER_T_HH

// ROOT includes
#include "Rtypes.h"

// EventNtuple includes
#include "EventNtuple/inc/TrigInfo.hh"

namespace Run1BEvtAna {
  struct Trigger_t  {
    const rooutil::Trigger* trig_;


    //-------------------------------------------------
    // Accessors
    bool Fired(int index) {
      if(!trig_) return false;
      return trig_->Fired(index);
    }

    bool Fired(const std::string& name) {
      if(!trig_) return false;
      return trig_->Fired(name);
    }

    //-------------------------------------------------
    // Additional functions

    bool FiredAPR() {
      if(!trig_) return false;
      for(const auto& itr : trig_->NameToIndexMap()) {
        const std::string& name = itr.first;
        bool found = false;
        found |= name.find("apr_TrkDe_75") != std::string::npos;
        found |= name.find("apr_TrkDe_80") != std::string::npos;
        if(found) {
          if(trig_->Fired(name)) return true;
        }
      }
      return false;
    }

    bool FiredCPR() {
      if(!trig_) return false;
      for(const auto& itr : trig_->NameToIndexMap()) {
        const std::string& name = itr.first;
        bool found = false;
        found |= name.find("cpr_TrkDe_75") != std::string::npos;
        found |= name.find("cpr_TrkDe_80") != std::string::npos;
        if(found) {
          if(trig_->Fired(name)) return true;
        }
      }
      return false;
    }

    void Reset() {
      trig_ = nullptr;
    }

    Trigger_t() { Reset(); }
  };
}
#endif
