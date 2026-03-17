#ifndef Run1BAnaUtils_hh
#define Run1BAnaUtils_hh
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

    //--------------------------------------------------------------------------------------
    static double closureApprox(const double energy, const double kmax = 90.1) {
      if(energy < 0. || energy > kmax || kmax == 0.) return 0.;
      const double x = energy/kmax;
      const double p = 20.*(1. - 2.*x + 2.*x*x)*x*(1.-x)*(1.-x);
      return p;
    }

    //--------------------------------------------------------------------------------------
    static CLHEP::HepLorentzVector lineAtCluster(const CaloCluster* cl, const KalSeed* seed,
                                                 GeomHandle<Calorimeter>& cal, const int debug_level = 0) {
      if(!cl || !seed) return CLHEP::HepLorentzVector(0.,0.,0., 0.);
      const int disk_id = cl->diskID();
      if(disk_id != 0 && disk_id != 1) return CLHEP::HepLorentzVector(0.,0.,0.,0.); // bad cluster
      if(seed->intersections().empty()) return CLHEP::HepLorentzVector(0.,0.,0.,0.); // bad line
      // take any intersection -- it's a line
      const auto line = seed->intersections().front().kinematicLine();
      const auto pos0 = line.pos0();
      const auto t0   = line.t0();
      const auto dir  = line.direction();
      if(std::abs(pos0.x()) > 1.e10) return CLHEP::HepLorentzVector(0.,0.,0., 0.); // bad line

      const auto pos_cl_trk = cal->geomUtil().mu2eToTracker(cal->geomUtil().diskToMu2e(disk_id, cl->cog3Vector()));
      const double dx = pos_cl_trk.x() - pos0.x();
      const double dy = pos_cl_trk.y() - pos0.y();
      const double dz = pos_cl_trk.z() - pos0.z();
      const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);

      // Get the position at the calo disk
      auto pos = pos0 + (dz/dir.z())*dir;
      CLHEP::Hep3Vector pos_trk(pos.x(), pos.y(), pos.z());
      const auto pos_cal = cal->geomUtil().mu2eToDisk(disk_id, cal->geomUtil().trackerToMu2e(pos_trk));

      // Get the time at the calo disk
      const double speed = line.speed();
      const double t_cal = t0 + ((pos0.z() < pos_cl_trk.z()) ? 1. : -1.) * distance / speed; // assume going ST -> calo

      // Validate the evaluation
      const double val_dz = pos_cal.z() - cl->cog3Vector().z();
      const bool val_issue =  std::abs(val_dz) > 1.;
      if(debug_level > 0 || val_issue) {
        std::cout << "[Run1BAnaUtils::" << __func__ << "]";
        if(val_issue) std::cout << " Problem with line matching!";
        std::cout << std::endl
                  << " Line pos0 = " << pos0 << std::endl
                  << " dir = " << dir << pos0 << std::endl
                  << " pos_trk = " << pos_trk << pos0 << std::endl
                  << " pos_cal = " << pos_cal << pos0 << std::endl
                  << " cluster = " << cl->cog3Vector() << pos0 << std::endl
                  << " t_cal = " << t_cal << pos0 << std::endl
                  << " t_cluster = " << cl->time() << pos0 << std::endl;
      }
      return CLHEP::HepLorentzVector(pos_cal.x(), pos_cal.y(), pos_cal.z(), t_cal);
    }

    //--------------------------------------------------------------------------------------
    static CLHEP::HepLorentzVector lineSeedAtCluster(const CaloCluster* cl, const CosmicTrackSeed* seed,
                                                     GeomHandle<Calorimeter>& cal, const int debug_level = 0) {
      if(!cl || !seed) return CLHEP::HepLorentzVector(0.,0.,0., 0.);
      const int disk_id = cl->diskID();
      if(disk_id != 0 && disk_id != 1) return CLHEP::HepLorentzVector(0.,0.,0.,0.); // bad cluster
      // retrieve the seed intercept and direction
      const auto seed_int = seed->track().Pos0();
      const auto seed_dir = seed->track().Dir().unit();
      const auto t0 = seed->t0().t0();

      const auto pos_cl_trk = cal->geomUtil().mu2eToTracker(cal->geomUtil().diskToMu2e(disk_id, cl->cog3Vector()));
      const double dx = pos_cl_trk.x() - seed_int.x();
      const double dy = pos_cl_trk.y() - seed_int.y();
      const double dz = pos_cl_trk.z() - seed_int.z();
      const double distance = std::sqrt(dx*dx + dy*dy + dz*dz);

      // Get the position at the calo disk

      auto pos_at_cal = seed_int + (dz/seed_dir.z())*seed_dir;
      CLHEP::Hep3Vector pos_trk(pos_at_cal.x(), pos_at_cal.y(), pos_at_cal.z());
      const auto pos_cal = cal->geomUtil().mu2eToDisk(disk_id, cal->geomUtil().trackerToMu2e(pos_trk));

      // Get the time at the calo disk
      const double speed = CLHEP::c_light; // assume c, good for electrons
      const double t_cal = t0 + distance / speed; // assume going ST -> calo

      if(debug_level > 0) {
        std::cout << "[Run1BAnaUtils::" << __func__ << "] " << std::endl
                  << " Seed int = " << seed_int << std::endl
                  << " dir = " << seed_dir << std::endl
                  << " pos_trk = " << pos_trk << std::endl
                  << " pos_cal = " << pos_cal << std::endl
                  << " cluster = " << cl->cog3Vector() << std::endl
                  << " t_cal = " << t_cal << std::endl
                  << " t_cluster = " << cl->time() << std::endl;
      }

      return CLHEP::HepLorentzVector(pos_cal.x(), pos_cal.y(), pos_cal.z(), t_cal);
    }

    //--------------------------------------------------------------------------------------
    // Printing utilities
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
        const auto& sim = sim_pair.second;
        out << std::format("{:8} {:8} {:8} {:8} {:10f} {:10f} {:10f} {:10f} {:10f} {:10f} {:10f} {:10f} {:20s}\n",
                           index, int(sim_pair.first.asUint()), (sim.hasParent()) ? sim.parent().key() : -1, int(sim.pdgId()),
                           sim.startMomentum().x(), sim.startMomentum().y(), sim.startMomentum().z(), sim.startMomentum().t(),
                           sim.startPosition().x(), sim.startPosition().y(), sim.startPosition().z(), sim.startGlobalTime(),
                           sim.creationCode().name().c_str());
        ++index;
      }
    }

    //--------------------------------------------------------------------------------------
    // Rough identification of sim particle type based on lineage
    enum {
      kUnknown = -1,
      kMuStop = 0,
      kMuStop_before_tracker = 1,
      kMuStop_after_tracker = 2,
      kMuStop_before_target = 3,
      kPrimary = 4,
      kNeutral = 5,
      kEleBeam = 6
    };
    static int getSimType(const SimParticle* sim) {
      if(!sim) return -1;
      int type = kUnknown;
      auto parent = sim->parent();
      while(parent.isNonnull()) {
        if(parent->pdgId() == 13) {
          if(sim->creationCode() == ProcessCode::mu2eMuonCaptureAtRest
             || sim->creationCode() == ProcessCode::mu2eMuonDecayAtRest) {
            type = kMuStop; // from a target muon stop
          } else if(parent->endPosition().z() < 4000.) { // before the target
            type = kMuStop_before_tracker;
          } else if(parent->endPosition().z() < 10000.) { // before the tracker
            type = kMuStop_before_tracker;
          } else if(parent->endPosition().z() > 10000.) { // after the tracker
            type = kMuStop_after_tracker;
          }
        } else if(!parent->hasParent()) { // last in the chain
          if(parent->creationCode() == ProcessCode::mu2ePrimary) {
            type = kPrimary;
          } else if(parent->creationCode() == ProcessCode::neutronInelastic) {
            type = kNeutral;
          } else if(std::abs(sim->pdgId()) != 13) {
            type = kEleBeam; // likely an electron from the beam, but could be something else
          }
        }
        if(type >= 0) break;
        sim = &*parent;
        parent = parent->parent();
      }
      return type;
    }
  };

} // namespace mu2e
#endif /* Run1BAnaUtils_hh */
