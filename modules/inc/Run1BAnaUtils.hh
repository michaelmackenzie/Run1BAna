//
//  Useful tools
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

// c++
#include <string>
#include <vector>
#include <iostream>

using namespace CLHEP;

namespace mu2e
{
  class Run1BAnaUtils
  {
  public:

    static void printSimCollection(const SimParticleCollection* sim_col, std::ostream& out = std::cout) {
      if(!sim_col) return;
      out << "Printing sim particle collection of size " << sim_col->size() << ":\n"
          << std::format("{:8s} {:8s} {:8s} {:8s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} {:10s} {:20s}\n",
                         "Index", "ID", "Parent", "Pdg ID",
                         "Start px", "Start py", "Start pz", "Start E",
                         "Start x", "Start y", "Start z", "Start t",
                         "Process code");
      size_t index = 0;
      for(const auto& sim_pair : *(sim_col)) {
        const auto& sim = *(sim_pair.second);
        out << "Printing sim particle collection of size " << sim_col->size() << ":\n"
            << std::format("{:8} {:8} {:8} {:8} {:10f} {:10f} {:10f} {:10f} {:10f} {:10f} {:10f} {:10f} {:20s}\n",
                           index, int(sim_pair.first.asUint()), (sim.hasParent()) ? sim.parent.key() : -1, int(sim.pdgId()),
                           sim.startMomentum().x(), sim.startMomentum().y(), sim.startMomentum().z(), sim.startMomentum().t(),
                           sim.startPosition().x(), sim.startPosition().y(), sim.startPosition().z(), sim.startPosition().t(),
                           sim.creationCode().name().c_str());
        ++index;
      }
    }

  };
}
