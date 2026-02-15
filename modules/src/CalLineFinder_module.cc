///////////////////////////////////////////////////////////////////////////////
// Calorimeter-driven straight line finding
// Pattern recognition stage
// Michael MacKenzie (2026)
// Based on CalHelixFinder, P. Murat and G. Pezzullo
///////////////////////////////////////////////////////////////////////////////

// framework
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// Offline
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrkReco/inc/TrkFaceData.hh"

#include "Offline/DataProducts/inc/Helicity.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixVal.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"

#include "Offline/Mu2eUtilities/inc/HelixTool.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"

// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BField/BField.hh"
#include "BTrk/BField/BFieldFixed.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BaBar/BbrStringUtils.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkErrCode.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
#include "BTrk/TrkBase/TrkParticle.hh"
#include "BTrk/TrkBase/TrkPoca.hh"

// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/algorithm/string.hpp>

//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"

// ROOT
#include "TROOT.h"
#include "TFolder.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TApplication.h"
#include "TFolder.h"
#include "TVector2.h"
#include "TSystem.h"
#include "TInterpreter.h"

// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <float.h>
#include <vector>
#include <set>
#include <map>

using namespace boost::accumulators;
using CLHEP::HepVector;
using CLHEP::Hep3Vector;

namespace mu2e {
  class Calorimeter;
  class Tracker;
  class ModuleHistToolBase;

  class CalLineFinder : public art::EDProducer {
  protected:

    //-----------------------------------------------------------------------------
    // Main module configuration parameters
    //-----------------------------------------------------------------------------
    struct Config
    {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>                           diag_level             {Name("DiagLevel"),                  Comment("Diagnostic output level"),0 };
      fhicl::Atom<art::InputTag>                 hit_coll_tag           {Name("ComboHitCollectionLabel"),    Comment(" Combo Hit Collection Label") };
      fhicl::Atom<std::string>                   time_cluster_coll_tag  {Name("TimeClusterCollectionLabel"), Comment("TimeCluster Collection Label") };
      fhicl::Atom<std::string>                   calo_cluster_coll_tag  {Name("CaloClusterCollectionLabel"), Comment("CaloCluster Collection Label") };
      fhicl::Atom<int>                           min_tc_hits            {Name("MinTimeClusterHits"),         Comment("Min NHits in TimeCluster") };
      fhicl::Atom<float>                         min_calo_cluster_energy{Name("MinCaloClusterEnergy"),       Comment("Min Calo Cluster Energy") };
      fhicl::Atom<float>                         max_edep_avg           {Name("MaxEDepAvg"),                 Comment("Max Avg EDep") };
      fhicl::Atom<int>                           min_line_hits          {Name("MinLineHits"),               Comment("Min NHits in Line") };
      fhicl::Atom<float>                         hit_time_sigma_thresh  {Name("HitTimeSigmaThresh"),         Comment("Time consistency threshold for hits to be added to the line") };
      fhicl::Atom<float>                         hit_xy_dist_thresh     {Name("HitXYDistThresh"),            Comment("Spatial consistency threshold for hits to be added to the line") };
      fhicl::Atom<std::string>                   fit_direction          {Name("FitDirection"),               Comment("Fit Direction in Search (\"downstream\" or \"upstream\")") };
      fhicl::Atom<bool>                          single_search          {Name("SingleSearch"),               Comment("Only search for one line per time cluster") };
    };

    //-----------------------------------------------------------------------------
    // Inputs
    //-----------------------------------------------------------------------------
    art::InputTag                         hit_tag_;                // input hit collection label
    art::InputTag                         time_cluster_tag_;       // input time cluster collection label
    art::InputTag                         calo_cluster_tag_;       // input calo cluster collection label

    int                                   min_tc_hits_;            // N(hits) in the time cluster
    float                                 min_calo_cluster_energy_;    // min energy of the associated calo cluster
    float                                 max_edep_avg_;           // max avg hit energy deposition
    size_t                                min_line_hits_;          // min number of hits in a line
    float                                 hit_time_sigma_thresh_;  // time consistency threshold for hits to be added to the line
    float                                 hit_xy_dist_thresh_;     // spatial consistency threshold for hits to be added to
    TrkFitDirection                       fit_dir_;                // fit direction in search
    bool                                  single_search_;          // only search for one line per time cluster (for diagnostics)
    int                                   diag_level_;             // diagnostic output
    //-----------------------------------------------------------------------------
    // Data
    //-----------------------------------------------------------------------------
    const art::Event*                     event_;

    const ComboHitCollection*             combo_hit_col_;           // input combo hit collection
    const TimeClusterCollection*          time_cluster_col_;        // input time cluster collection
    const CaloClusterCollection*          calo_cluster_col_;        // input calo cluster collection
    art::Handle<ComboHitCollection>       combo_hit_col_handle_;    // handle for input combo hit collection
    art::Handle<TimeClusterCollection>    time_cluster_col_handle_; // handle for input time cluster collection

    const Tracker*                        tracker_     ;            // tracker geometry
    const Calorimeter*                    calorimeter_ ;            // calorimeter geometry
    CLHEP::Hep3Vector                     target_pos_;              // position of the target center for seeding
    double                                calo_d0_offset_;          // z offset of the calorimeter disk 0 from the tracker system
    double                                calo_d1_offset_;          // z offset of the calorimeter disk 1 from the tracker system
    std::unordered_set<size_t>            hits_used_in_lines_;       // set of hit indices that have already been used in lines (for single line search)
  public:
    explicit CalLineFinder(const art::EDProducer::Table<Config>& config);
    virtual ~CalLineFinder();

    virtual void beginJob();
    virtual void beginRun(art::Run&   run   );
    virtual void produce (art::Event& event );
    virtual void endJob();
    //-----------------------------------------------------------------------------
    // helper functions
    //-----------------------------------------------------------------------------
    bool findData             (const art::Event& e);
    bool isGoodTimeCluster    (const TimeCluster& tc);
    int  goodHitsTimeCluster  (const TimeCluster& TimeCluster);
    bool findLineInTimeCluster(const size_t i_tc, std::vector<CosmicTrackSeed>& seeds);
    void fitLine(const CLHEP::Hep3Vector& cluster_pos, CLHEP::Hep3Vector& seed_dir, CLHEP::Hep3Vector& seed_int, std::vector<size_t>& hit_indices);
  };

  //-----------------------------------------------------------------------------
  // module constructor, parameter defaults are defiend in CalPatRec/fcl/prolog.fcl
  //-----------------------------------------------------------------------------
  CalLineFinder::CalLineFinder(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config}
    , hit_tag_(config().hit_coll_tag())
    , time_cluster_tag_(config().time_cluster_coll_tag())
    , calo_cluster_tag_(config().calo_cluster_coll_tag())
    , min_tc_hits_(config().min_tc_hits())
    , min_calo_cluster_energy_(config().min_calo_cluster_energy())
    , max_edep_avg_(config().max_edep_avg())
    , min_line_hits_(config().min_line_hits())
    , hit_time_sigma_thresh_(config().hit_time_sigma_thresh())
    , hit_xy_dist_thresh_(config().hit_xy_dist_thresh())
    , fit_dir_(config().fit_direction())
    , diag_level_(config().diag_level())
   {
    // declare the data products
    consumes<ComboHitCollection>(hit_tag_);
    consumes<TimeClusterCollection>(time_cluster_tag_);
    consumes<CaloClusterCollection>(calo_cluster_tag_);
    produces<CosmicTrackSeedCollection>();
  }

  //-----------------------------------------------------------------------------
  // destructor
  //-----------------------------------------------------------------------------
  CalLineFinder::~CalLineFinder() {
  }

  //-----------------------------------------------------------------------------
  void CalLineFinder::beginJob() {
    // art::ServiceHandle<art::TFileService> tfs;
  }

  //-----------------------------------------------------------------------------
  void CalLineFinder::beginRun(art::Run& ) {
    mu2e::GeomHandle<mu2e::Tracker> th;
    tracker_ = th.get();

    mu2e::GeomHandle<mu2e::Calorimeter> ch;
    calorimeter_ = ch.get();

    // get the offset between the calo disks and the tracker system
    calo_d0_offset_ = calorimeter_->geomUtil().mu2eToTracker(calorimeter_->geomUtil().diskToMu2e(0, CLHEP::Hep3Vector(0., 0., 0.))).z();
    calo_d1_offset_ = calorimeter_->geomUtil().mu2eToTracker(calorimeter_->geomUtil().diskToMu2e(1, CLHEP::Hep3Vector(0., 0., 0.))).z();
  }

  //-----------------------------------------------------------------------------
  // find the input data objects
  //-----------------------------------------------------------------------------
  bool CalLineFinder::findData(const art::Event& evt) {

    evt.getByLabel(hit_tag_, combo_hit_col_handle_);
    if(!combo_hit_col_handle_.isValid()) {
      printf("[CalLineFinder::%s] ERROR: ComboHit collection with label \"%s\" not found! RETURN\n", __func__, hit_tag_.encode().c_str());
      combo_hit_col_ = nullptr;
      return false;
    } else {
      combo_hit_col_ = combo_hit_col_handle_.product();
    }

    auto calo_cluster_handle = evt.getValidHandle<CaloClusterCollection>(calo_cluster_tag_);
    calo_cluster_col_ = calo_cluster_handle.product();

    evt.getByLabel(time_cluster_tag_, time_cluster_col_handle_);
    if(time_cluster_col_handle_.isValid()) {
      time_cluster_col_ = time_cluster_col_handle_.product();
    } else {
      printf("[CalLineFinder::%s] ERROR: TimeCluster collection with label \"%s\" not found! RETURN\n", __func__, time_cluster_tag_.encode().c_str());
      time_cluster_col_ = nullptr;
      return false;
    }

    return true;
  }

  //-----------------------------------------------------------------------------
  // Find lines in a time cluster and add them to the output seed collection
  //-----------------------------------------------------------------------------
  bool CalLineFinder::findLineInTimeCluster(const size_t tc_index, std::vector<CosmicTrackSeed>& seeds) {
    const auto& tc = time_cluster_col_->at(tc_index);
    // Get the associated calo cluster and use its position/time as a seed for line finding
    const auto& cl = *tc.caloCluster();
    const float cl_time = cl.time();
    const bool downstream = fit_dir_ == TrkFitDirection::FitDirection::downstream;
    CosmicTrackSeed seed;

    // // get collection at straw level
    // std::vector<StrawHitIndex> shiv;
    // auto hit_col_ptr = combo_hit_col_->fillStrawHitIndices(tc.hits(), shiv, StrawIdMask::uniquestraw);
    // auto combo_hits = *hit_col_ptr;
    // seed._straw_chits.setAsSubset(combo_hit_col_handle_, StrawIdMask::uniquestraw);


    CLHEP::Hep3Vector cl_pos = cl.cog3Vector();
    cl_pos.setZ(cl_pos.z() + (cl.diskID() == 0 ? calo_d0_offset_ : calo_d1_offset_)); // shift the cluster position to the tracker system

    CLHEP::Hep3Vector seed_dir = (downstream) ? cl_pos - target_pos_ : target_pos_ - cl_pos;
    seed_dir.setX(seed_dir.x()/2.); // assume some slope, but not too steep
    seed_dir.setY(seed_dir.y()/2.);
    const double seed_dir_mag = seed_dir.mag();
    if(seed_dir_mag <= 0.) return false; // can't define a seed direction, so return false{
    seed_dir *= 1./seed_dir_mag; // normalize the seed direction
    if(diag_level_ > 1) {
      printf("[CalLineFinder::%s] TimeCluster %zu: CaloCluster time = %.2f, position = (%.1f, %.1f, %.1f), seed direction = (%.2f, %.2f, %.2f)\n",
             __func__, tc_index, cl_time, cl_pos.x(), cl_pos.y(), cl_pos.z(), seed_dir.x(), seed_dir.y(), seed_dir.z());
    }

    // LsqSums2 line_fitter; // fitter for finding lines
    // Look for hits in the time cluster that are consistent with this seed position and direction, and add them to the seed
    std::vector<size_t> hit_indices;
    hit_indices.reserve(tc.hits().size());
    for(size_t i_hit = 0; i_hit < tc.hits().size(); ++i_hit) {
      const size_t hit_index = tc.hits().at(i_hit);
      if(hits_used_in_lines_.count(hit_index) > 0) continue; // this hit has already been used in a line, so skip it
      const auto& hit = combo_hit_col_->at(hit_index);
      const CLHEP::Hep3Vector hit_pos(hit.pos().x(), hit.pos().y(), hit.pos().z());
      const double hit_time = hit.correctedTime();
      const double time_at_hit = cl_time + (hit_pos - cl_pos).dot(seed_dir) / CLHEP::c_light; // time at the hit position based on the seed direction and cluster time
      const double time_sigma = std::abs((hit_time - time_at_hit) / hit.timeRes()); // number of sigma the hit time is from the expected time based on the seed
      if(diag_level_ > 2) {
        printf(" Hit %zu: time = %.2f, expected time = %.2f, time sigma = %.2f\n",
               hit_index, hit_time, time_at_hit, time_sigma);
      }
      if(time_sigma > hit_time_sigma_thresh_) continue; // hit is not consistent with the seed, so skip it

      // next check if the hit is consistent in space
      const double dz = (seed_dir.z() > 0.) ? hit_pos.z() - cl_pos.z() : cl_pos.z() - hit_pos.z(); // z distance from the cluster to the hit along the seed direction
      const CLHEP::Hep3Vector pos_at_hit = cl_pos + dz * seed_dir; // position along the seed direction at the same z as the hit
      const double x_y_dist = (hit_pos - pos_at_hit).perp(); // distance in the x-y plane between the hit and the expected position based on the seed
      if(diag_level_ > 2) {
        printf("  Expected position at hit z: (%.1f, %.1f, %.1f), hit position: (%.1f, %.1f, %.1f), x-y distance = %.2f\n",
               pos_at_hit.x(), pos_at_hit.y(), pos_at_hit.z(), hit_pos.x(), hit_pos.y(), hit_pos.z(), x_y_dist);
      }
      if(x_y_dist > hit_xy_dist_thresh_) continue; // hit is not consistent with the seed
      hit_indices.push_back(hit_index);
    }

    if(diag_level_ > 1) {
      printf("[CalLineFinder::%s] TimeCluster %zu: Found %zu hits consistent with the seed\n", __func__, tc_index, hit_indices.size());
    }
    if(hit_indices.size() < min_line_hits_) return false; // not enough hits consistent with the seed, so return false

    // fit the line parameters based on the cluster position and the consistent hits
    CLHEP::Hep3Vector seed_int; // point of closest approach to the z axis (x0, y0, 0)
    fitLine(cl_pos, seed_dir, seed_int, hit_indices);

      // get pos and direction into Z alignment
    if (seed_dir.y() != 0){
      seed_dir /= -1*seed_dir.y();
      seed_int -= seed_dir*seed_int.y()/seed_dir.y();
    }

    if(diag_level_ > 0) {
      printf("[CalLineFinder::%s] TimeCluster %zu: Fitted line parameters: seed_int = (%.1f, %.1f, %.1f), seed_dir = (%.2f, %.2f, %.2f)\n",
             __func__, tc_index, seed_int.x(), seed_int.y(), seed_int.z(), seed_dir.x(), seed_dir.y(), seed_dir.z());
    }

    for(size_t hit_index : hit_indices) {
      hits_used_in_lines_.insert(hit_index); // add this hit index to the set of hits used in lines (for single line search)
      ComboHit hit(combo_hit_col_->at(hit_index)); // clone the combo hit
      seed._straw_chits.push_back(std::move(hit));
    }
    seed._straw_chits.setParent(combo_hit_col_->parent());

    // // double avg_t0 = 0.;
    // // int good_hits = 0;
    // // seed._straw_chits.setAsSubset(combo_hit_col_handle_,StrawIdMask::uniquestraw);
    // for (size_t k=0;k<shiv.size();k++){
    //   size_t kloc = shiv[k];
    //   // double traj_time = 0.;
    //   // double hit_t0 = shC[kloc].time() - shC[kloc].driftTime() - shC[kloc].propTime() - traj_time;
    //   // avg_t0 += hit_t0;
    //   // good_hits++;
    //   ComboHit combohit;
    //   combohit.init(combo_hit_col_[kloc],kloc);
    //   seed._straw_chits.push_back(std::move(combohit));
    // }

    // avg_t0 /= good_hits;

    seed._t0._t0 = cl_time; // set the seed t0 to the cluster time
    seed._t0._t0err = 1.; // dummy error on seed t0
    seed._timeCluster = art::Ptr<TimeCluster>(time_cluster_col_handle_, tc_index);
    seed._caloCluster = tc.caloCluster();
    seed._track.converged = true;
    seed._track.FitParams.T0 = seed._t0._t0;
    seed._track.FitParams.A0 = seed_int.x();
    seed._track.FitParams.B0 = seed_int.z();
    seed._track.FitParams.A1 = seed_dir.x();
    seed._track.FitParams.B1 = seed_dir.z();
    seed._track.MinuitParams.T0 = seed._t0._t0;
    seed._track.MinuitParams.A0 = seed_int.x();
    seed._track.MinuitParams.B0 = seed_int.z();
    seed._track.MinuitParams.A1 = seed_dir.x();
    seed._track.MinuitParams.B1 = seed_dir.z();
    XYZVectorF X(1,0,0);
    XYZVectorF Y(0,1,0);
    XYZVectorF Z(0,0,1);
    TrackAxes XYZ(X,Y,Z);
    seed._track.InitCoordSystem = XYZ;
    seed._track.FitCoordSystem = XYZ;
    XYZVectorF xyzint(seed_int.x(), seed_int.y(), seed_int.z());
    XYZVectorF xyzdir(seed_dir);
    TrackEquation XYZTrack(xyzint,xyzdir);
    seed._track.SetFitEquation(XYZTrack);
    seed._track.SetMinuitEquation(XYZTrack);
    seed._status.merge(TrkFitFlag::Straight);
    seed._status.merge(TrkFitFlag::hitsOK);
    seed._status.merge(TrkFitFlag::helixOK);
    seed._status.merge(TrkFitFlag::helixConverged);
    seed._track.MinuitParams.cov = std::vector<double>(15, 0);
    seeds.push_back(seed);
    return true; // for now, just add the seed based on the calo cluster position and return false to keep searching for more lines in this time cluster (if single_search_ is false)
  }

  void CalLineFinder::fitLine(const CLHEP::Hep3Vector& cluster_pos,
                              CLHEP::Hep3Vector& seed_dir,
                              CLHEP::Hep3Vector& seed_int,
                              std::vector<size_t>& hit_indices) {
    // Fix one end of the line to the cluster position and fit the dx/dz and dy/dz slopes based on the hits
    // x(z) = m_x * (z - z0) + x0 --> z(x) = (x - x0)/m_x + z0
    // y(z) = m_y * (z - z0) + y0 --> y(x) = m_y * ((x - x0)/m_x + z0 - z0) + y0 = (m_y/m_x) * (x - x0) + y0
    // where (x0, y0, z0) is the cluster position and m_x = dx/dz and m_y = dy/dz are the slopes to be fitted

    // if dz is 0, then we can't fit the line, so just return
    if(std::fabs(seed_dir.z()) < 0.01) return;

    // need at least 2 hits + calo hit
    if(hit_indices.size() < 2) return;

    double       m_x = seed_dir.x() / seed_dir.z(); // initial guess for dx/dz based on the seed direction
    double       m_y = seed_dir.y() / seed_dir.z(); // initial guess for dy/dz based on the seed direction
    const double z0  = cluster_pos.z();
    double       x0  = cluster_pos.x();
    double       y0  = cluster_pos.y();

    LsqSums2 x_z_fitter(z0, x0); // fitter for x vs z
    LsqSums2 y_z_fitter(z0, y0); // fitter for y vs z
    // add the cluster location with a small weight to the fit to help constrain it, but not too much that it dominates over the hits
    x_z_fitter.addPoint(z0, x0, 0.1);
    y_z_fitter.addPoint(z0, y0, 0.1);
    for(size_t hit_index : hit_indices) {
      const auto& hit = combo_hit_col_->at(hit_index);
      const double z = hit.pos().z();
      const double x = hit.pos().x();
      const double y = hit.pos().y();
      x_z_fitter.addPoint(z, x, 1.0);
      y_z_fitter.addPoint(z, y, 1.0);
    }
    if(x_z_fitter.qn() > 1) {
      m_x = x_z_fitter.dydx();
      x0  = x_z_fitter.y0();
    }
    if(y_z_fitter.qn() > 1) {
      m_y = y_z_fitter.dydx();
      y0  = y_z_fitter.y0();
    }
    seed_dir.setX(m_x * seed_dir.z());
    seed_dir.setY(m_y * seed_dir.z());
    seed_dir *= 1./seed_dir.mag(); // normalize the seed direction after fitting

    seed_int.setX(x0 - m_x * z0);
    seed_int.setY(y0 - m_y * z0);
    if(diag_level_ > 1) {
      printf("[CalLineFinder::%s] Fitted line parameters: seed_int = (%.1f, %.1f, %.1f), seed_dir = (%.2f, %.2f, %.2f)\n",
             __func__, seed_int.x(), seed_int.y(), seed_int.z(), seed_dir.x(), seed_dir.y(), seed_dir.z());
    }
  }

  //-----------------------------------------------------------------------------
  // event entry point
  //-----------------------------------------------------------------------------
  void CalLineFinder::produce(art::Event& event ) {

    // diagnostic info
    event_     = &event;

    // output collection
    std::unique_ptr<CosmicTrackSeedCollection> out_seeds(new CosmicTrackSeedCollection);

    //-----------------------------------------------------------------------------
    // find the data
    //-----------------------------------------------------------------------------

    const bool valid_data = findData(event);

    //-----------------------------------------------------------------------------
    // Search in each time cluster for line candidates
    //-----------------------------------------------------------------------------

    if(valid_data) {
      const size_t n_time_clusters = time_cluster_col_->size();
      for(size_t i_tc = 0; i_tc < n_time_clusters; ++i_tc) {
        const auto& tc = time_cluster_col_->at(i_tc);
        if(!isGoodTimeCluster(tc)) continue; // only search in good time clusters

        // Search for lines in this time cluster
        std::vector<CosmicTrackSeed> seeds;
        while(findLineInTimeCluster(i_tc, seeds) && !single_search_) { // keep searching until no more lines are found in this time cluster
          if(diag_level_ > 1) {
            printf("[CalLineFinder::%s] Found line candidate in time cluster %zu with %zu hits\n", __func__, i_tc, seeds.back()._straw_chits.size());
          }
        }
        if(diag_level_ > 0) {
          printf("[CalLineFinder::%s] Found %zu line candidates in time cluster %zu (n hits in time cluster = %zu)\n",
                 __func__, seeds.size(), i_tc, tc.nhits());
        }
        if(seeds.empty()) continue;
        out_seeds->insert(out_seeds->end(), seeds.begin(), seeds.end());
      }
    } else {
      printf("[CalLineFinder::%s] ERROR: Input data not found! RETURN\n", __func__);
    }

    //-----------------------------------------------------------------------------
    // put reconstructed seeds into the event record
    //-----------------------------------------------------------------------------

    if(diag_level_ > 0) {
      printf("[CalLineFinder::%s] Found %zu line candidates in total\n", __func__, out_seeds->size());
    }
    event.put(std::move(out_seeds));
  }

  int  CalLineFinder::goodHitsTimeCluster(const TimeCluster& TimeCluster){
    int   nhits         = TimeCluster.nhits();
    int   ngoodhits(0);
    for (int i=0; i<nhits; ++i){
      const int    index   = TimeCluster.hits().at(i);
      const auto&  hit     = combo_hit_col_->at(index);
      const auto&  flag    = hit.flag();
      const int    bkg_hit = flag.hasAnyProperty(StrawHitFlag::bkg);
      if (bkg_hit != 0) continue;
      ngoodhits += hit.nStrawHits();
    }

    return ngoodhits;
  }

  bool CalLineFinder::isGoodTimeCluster(const TimeCluster& tc) {
    if(!tc.hasCaloCluster()) return false; // must have an associated calorimeter cluster
    const auto& cl = *tc.caloCluster();
    if(cl.energyDep() < min_calo_cluster_energy_) return false; // must have at least min_calo_energy_ energy deposition
    if(goodHitsTimeCluster(tc) < min_tc_hits_) return false; // must have at least min_tc_hits_ hits
    return true;
  }
  //-----------------------------------------------------------------------------
  //
  //-----------------------------------------------------------------------------
  void CalLineFinder::endJob() {
  }

} // end namespace mu2e

using mu2e::CalLineFinder;
DEFINE_ART_MODULE(CalLineFinder)
