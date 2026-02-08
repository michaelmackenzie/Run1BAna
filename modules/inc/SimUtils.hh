//
//  Help associate sim particles to straw hits/digis based on MC truth
//  Michael MacKenzie, 2026
//

// Mu2e
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/MCDataProducts/inc/CaloClusterMC.hh"
#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include "Offline/MCDataProducts/inc/KalSeedMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"

// ROOT
#include "TH1.h"
#include "TString.h"

// c++
#include <string>
#include <vector>
#include <iostream>

using namespace CLHEP;

namespace mu2e
{
  class SimUtils
  {
  public:

    //--------------------------------------------------------------------------------------x
    enum SimType {
      kElectron = 0, kProton, kPhoton, kNeutron, kMuon, kPion,
      kIPA, kMuElectron, kMuCapProton, kMuCapNeutron, kMuCapPhoton,
      kIon, kOther, kNone, kUnknown, kLast
    };

    //--------------------------------------------------------------------------------------
    // For sim particle info
    struct Sim_t {
      unsigned id_ = 0;
      unsigned nhits_ = 0;
      float hit_start_z_ = 0.;
      float hit_start_p_ = 0.;
      float hit_start_t_ = 0.;
      float hit_end_z_ = 0.;
      float hit_end_p_ = 0.;
      float hit_end_t_ = 0.;
      float start_p_ = 0.;
      float start_pz_ = 0.;
      float start_z_ = 0.;
    };

    static void fillSimInfo(std::map<unsigned, Sim_t>& simInfo, const SimParticleCollection* simcol, const StrawDigiMCCollection* digcol);

    static const SimParticle* simByID(const SimParticleCollection* simcol, const unsigned id) {
      return simcol->getOrNull(cet::map_vector_key(id));
    }

    static bool checkIfParent(const SimParticleCollection* simcol, const SimParticle* sim_1, const SimParticle* sim_2) {
      if(!sim_1 || !sim_2) return false;
      if(sim_1->id() == sim_2->id()) return true; // self-parent-like

      // Check if sim 1 is a parent to sim 2
      art::Ptr<SimParticle> parent = sim_2->parent();
      while(parent.isNonnull()) {
        if(sim_1->id() == parent->id()) return true; // parent relationship found
        parent = parent->parent();
      }

      // Check if sim 2 is a parent to sim 1
      parent = sim_1->parent();
      while(parent.isNonnull()) {
        if(sim_2->id() == parent->id()) return true; // parent relationship found
        parent = parent->parent();
      }
      return false;
    }

    static bool checkIfParent(const SimParticleCollection* simcol, const unsigned id_1, const unsigned id_2) {
      return checkIfParent(simcol, simByID(simcol, id_1), simByID(simcol, id_2));
    }

    static bool isIPA(const SimParticle* sim) {
      if(!sim) return false;
      const SimParticle* parent = (sim->hasParent()) ? &(*sim->parent()) : nullptr;
      if(parent && std::abs(parent->pdgId()) == 13) {
        const float x(parent->endPosition().x()+3904.),
          y(parent->endPosition().y()), z(parent->endPosition().z());
        const float r(std::sqrt(x*x+y*y));
        return std::fabs(r - 300.) < 10. && z > 6500. && z < 10000.; //  FIXME: Get z limits
      }
      return false;
    }

    static bool isIPAOrigin(const SimParticle* sim) {
      if(!sim) return false;
      // continue up the lineage checking if it fails
      return isIPA(sim) || (sim->parent().isNonnull() && isIPAOrigin(&(*sim->parent())));
    }

    static SimType getSimType(const SimParticle* sim) {
      if(!sim) return SimType::kUnknown;
      if(isIPA(sim)) return SimType::kIPA;
      const int pdg = sim->pdgId();
      const bool mu_daughter = (sim->hasParent() && sim->parent().isAvailable()
                                && std::abs(sim->parent()->pdgId()) == 13)
        || sim->creationCode() == ProcessCode::mu2eMuonCaptureAtRest
        || sim->creationCode() == ProcessCode::mu2eMuonDecayAtRest;
      if(mu_daughter) {
        switch(pdg) {
        case 11: case -11:   return SimType::kMuElectron;
        case 2212:           return SimType::kMuCapProton;
        case 2112:           return SimType::kMuCapNeutron;
        case 22:             return SimType::kMuCapPhoton;
        }
      }
      switch(pdg) {
      case 11: case -11:   return SimType::kElectron;
      case 13: case -13:   return SimType::kMuon;
      case 211: case -211: return SimType::kPion; // pi+/-
      case 111:            return SimType::kPion; // pi0
      case 2212:           return SimType::kProton;
      case 2112:           return SimType::kNeutron;
      case 22:             return SimType::kPhoton;
      }
      if(std::abs(pdg) > 10000) return SimType::kIon;
      std::cout << __func__ << ": Unknown PDG: " << pdg << std::endl;
      return SimType::kOther;
    }
    static TString getSimTypeName(SimType type) {
      switch(type) {
      case SimType::kElectron: return "Electron";
      case SimType::kProton  : return "Proton";
      case SimType::kPhoton  : return "Photon";
      case SimType::kNeutron : return "Neutron";
      case SimType::kMuElectron   : return "DIO";
      case SimType::kMuCapProton  : return "#mu-cap p^{+}";
      case SimType::kMuCapPhoton  : return "#mu-cap #gamma";
      case SimType::kMuCapNeutron : return "#mu-cap n^{0}";
      case SimType::kMuon    : return "Muon";
      case SimType::kPion    : return "Pion";
      case SimType::kIPA     : return "IPA DIO";
      case SimType::kIon     : return "Ion";
      case SimType::kOther   : return "Other";
      case SimType::kNone    : return "None";
      default: return "Unknown";
      }
    }
    static TString getSimTypeName(const SimParticle* sim) {
      const auto type = getSimType(sim);
      return getSimTypeName(type);
    }
    TString getSimTypeName(const art::Ptr<SimParticle> simptr) {
      const SimParticle* sim = (simptr.isNonnull()) ? &(*simptr) : nullptr;
      const auto type = getSimType(sim);
      return getSimTypeName(type);
    }

    static SimType getSimOriginType(const SimParticle* sim, const int depth = 0, int debugLevel = 0);

  };


  //--------------------------------------------------------------------------------------
  void SimUtils::fillSimInfo(std::map<unsigned, Sim_t>& simInfo,
                             const SimParticleCollection* simcol, const StrawDigiMCCollection* digcol) {
    simInfo.clear(); // clear the info
    if(!simcol || !digcol) return;

    // Count the number of hits for each sim particle
    struct info_t {
      unsigned nhits_ = 0;
      float hit_start_z_ = 1.e10;
      float hit_start_p_ = 0.;
      float hit_start_t_ = 0.;
      float hit_end_z_ = -1.e10;
      float hit_end_p_ = 0.;
      float hit_end_t_ = 0.;
    };
    std::map<unsigned int,unsigned> pmap;
    std::map<unsigned int, info_t> info_map;
    for(auto const& digi : *digcol) {
      // do not inspect truth info for digis not produced in simulation
      if (digi.containsSimulation()){
        // look at the early end
        StrawEnd fend = digi.earlyEnd();
        auto const& step =  digi.strawGasStep(fend);
        art::Ptr<SimParticle> const& sp = step->simParticle();
        const unsigned id = sp->id().asInt();
        ++pmap[id];

        const float z = step->position().z();
        const float p = std::sqrt(step->momentum().mag2());
        const float t = step->time();
        // if(debugLevel > 5) printf(" ID = %2i, z = %7.1f, p = %6.2f, t = %6.1f\n",
        //                            (int) id, z, p, t);
        if(info_map.find(id) == info_map.end()) info_map[id] = info_t();
        info_t& info = info_map[id];
        if(info.hit_start_z_ > z) {
          info.hit_start_z_ = z;
          info.hit_start_p_ = p;
          info.hit_start_t_ = t;
        }
        if(info.hit_end_z_ < z) {
          info.hit_end_z_ = z;
          info.hit_end_p_ = p;
          info.hit_end_t_ = t;
        }
      }
    }

    // Store the information for the sim particles of interest
    for(auto sim : *simcol) {
      Sim_t sim_t;
      sim_t.id_ = sim.second.id().asInt();
      auto mapfnd = pmap.find(sim_t.id_);
      sim_t.nhits_ = (mapfnd == pmap.end()) ? 0 : mapfnd->second;
      if(mapfnd != pmap.end()) {
        const info_t& info = info_map[sim_t.id_];
        sim_t.hit_start_z_ = info.hit_start_z_;
        sim_t.hit_start_p_ = info.hit_start_p_;
        sim_t.hit_start_t_ = info.hit_start_t_;
        sim_t.hit_end_z_   = info.hit_end_z_  ;
        sim_t.hit_end_p_   = info.hit_end_p_  ;
        sim_t.hit_end_t_   = info.hit_end_t_  ;
      }
      simInfo[sim_t.id_] = sim_t;
    }
  }

  //--------------------------------------------------------------------------------------
  // Determine the source process for this sim
  SimUtils::SimType SimUtils::getSimOriginType(const SimParticle* sim, int depth, int debugLevel) {
    if(!sim) return SimType::kUnknown;

    const SimParticle* current = sim;
    while(current) {
      ++depth;
      art::Ptr<SimParticle> parent = current->parent();
      if(debugLevel > 3) std::cout << std::string(depth, ' ') << " level " << depth << " pdg " << current->pdgId()
                                   << " creation " << current->creationCode()
                                   << " stopping " << current->stoppingCode()
                                   << " start x " << current->startPosition()
                                   << " and p " << current->startMomentum().vect().mag()
                                   << " has parent " << current->hasParent()
                                   << std::endl;

      // Check if the parent originates from the IPA
      if(parent.isAvailable() && isIPA(&(*parent))) {
        if(debugLevel > 3) std::cout << " --> From IPA muon\n";
        return SimType::kIPA;
      }

      // Check if the parent is a muon from the stopping target
      if(parent.isAvailable() && std::abs(parent->pdgId()) == 13) { // muon parent
        if(parent->endPosition().z() > 4500. &&
           parent->endPosition().z() < 7000.) { // target muon stop
          auto type = getSimType(current);
          if(debugLevel > 3) std::cout << " --> From ST muon, Type "
                                       << getSimTypeName(type).Data() << std::endl;
          return type;
        }
      }

      // Check if a standard muon stop creation code process
      if(current->creationCode() == ProcessCode::mu2eMuonCaptureAtRest
         || current->creationCode() == ProcessCode::mu2eMuonDecayAtRest) {
          auto type = getSimType(current);
          if(debugLevel > 3) std::cout << " --> From a muon, Type "
                                       << getSimTypeName(type).Data() << std::endl;
          return type;
      }


      // next level up, if available
      if(!parent.isAvailable() || parent.isNull()) break;
      current = &(*parent);
    }

    if(current) {
      auto type = getSimType(current);
      if(debugLevel > 3) std::cout << " --> Type " << getSimTypeName(type).Data() << std::endl;
      return type;
    }
    auto type = getSimType(sim);
    if(debugLevel > 3) std::cout << " --> No history, Type " << getSimTypeName(type).Data() << std::endl;
    return type;
  }

}
