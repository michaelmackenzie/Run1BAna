//
// Sim particle information
// Michael MacKenzie (2025)

#ifndef RUN1BEVTANA_SIMPARTICLE_T_HH
#define RUN1BEVTANA_SIMPARTICLE_T_HH

// ROOT includes
#include "Rtypes.h"
#include "Math/LorentzVector.h"

// Event ntuple includes
#include "EventNtuple/inc/SimInfo.hh"

// Local includes
#include "Run1BAna/evtana/inc/GlobalConstants.h"

namespace Run1BEvtAna {

  struct SimParticle_t {
    ROOT::Math::XYZTVectorF mom_start_;
    ROOT::Math::XYZTVectorF mom_end_  ;
    ROOT::Math::XYZTVectorF pos_start_;
    ROOT::Math::XYZTVectorF pos_end_  ;
    Int_t id_        ;
    Int_t pdg_       ;
    Int_t start_code_;
    Int_t end_code_  ;
    Int_t nhits_     ;
    Int_t ntrkhits_  ;
    Int_t mcrel_     ;

    void Initialize(mu2e::SimInfo* info) {
      if(!info) { Reset(); return; }
      const int pdg = info->pdg;
      // auto const& ptable = GlobalConstantsHandle<ParticleDataList>();
      // const double mass_ = ptable->particle(fpart_).mass();
      const double mass((std::abs(pdg) > 10000) ? 0. : ParticleMass(pdg));
      const float start_e = std::sqrt(std::pow(info->   mom.x(), 2) + std::pow(info->   mom.y(),2) + std::pow(info->   mom.z(), 2) + mass*mass);
      const float end_e   = std::sqrt(std::pow(info->endmom.x(), 2) + std::pow(info->endmom.y(),2) + std::pow(info->endmom.z(), 2) + mass*mass);
      mom_start_.SetPxPyPzE(info->mom.x(), info->mom.y(), info->mom.z(), start_e);
      mom_end_  .SetPxPyPzE(info->endmom.x(), info->endmom.y(), info->endmom.z(), end_e);
      pos_start_.SetXYZT(info->pos.x(), info->pos.y(), info->pos.z(), info->time);
      pos_end_  .SetXYZT(info->endpos.x(), info->endpos.y(), info->endpos.z(), info->time + 1e6); // FIXME: No end time available
      id_         = info->id;
      pdg_        = info->pdg;
      start_code_ = info->startCode;
      end_code_   = info->stopCode;
      nhits_      = info->nhits;
      ntrkhits_   = info->nactive;
      mcrel_      = info->prirel.relationship();
    }

    void Reset() {
      mom_start_.SetXYZT(0., 0., 0., 0.);
      mom_end_  .SetXYZT(0., 0., 0., 0.);
      pos_start_.SetXYZT(0., 0., 0., 0.);
      pos_end_  .SetXYZT(0., 0., 0., 0.);
      id_         = 0;
      pdg_        = 0;
      start_code_ = 0;
      end_code_   = 0;
      nhits_      = 0;
      ntrkhits_   = 0;
      mcrel_      = -1000;
    }
    SimParticle_t() { Reset(); }
  };
}
#endif
