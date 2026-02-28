//
//  Study Run 1B info
//  Michael MacKenzie, 2026
//
// framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/System/TriggerNamesService.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/TriggerResults.h"

// Mu2e
#include "Offline/CaloCluster/inc/ClusterUtils.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/MCDataProducts/inc/CaloClusterMC.hh"
#include "Offline/MCDataProducts/inc/CaloShowerSim.hh"
#include "Offline/MCDataProducts/inc/GenEventCount.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/Mu2eUtilities/inc/StopWatch.hh"
#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/CaloGeomUtil.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

// MC truth associations
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"

// ROOT
#include "TH1.h"
#include "TH2.h"
#include "TString.h"

// c++
#include <string>
#include <vector>
#include <iostream>

// local
#include "Run1BAna/modules/inc/Structs.hh"
#include "Run1BAna/modules/inc/SimUtils.hh"
#include "Run1BAna/modules/inc/Run1BAnaUtils.hh"

using namespace CLHEP;
using EventHist_t = Run1BAnaStructs::EventHist_t;
using ClusterHist_t = Run1BAnaStructs::ClusterHist_t;
using LineHist_t = Run1BAnaStructs::LineHist_t;
using CosmicSeedHist_t = Run1BAnaStructs::CosmicSeedHist_t;
using TimeClusterHist_t = Run1BAnaStructs::TimeClusterHist_t;
using SimHist_t = Run1BAnaStructs::SimHist_t;
// using ChargedTrackHist_t = Run1BAnaStructs::ChargedTrackHist_t;

using EventPar_t = Run1BAnaStructs::EventPar_t;
using ClusterPar_t = Run1BAnaStructs::ClusterPar_t;
using LinePar_t = Run1BAnaStructs::LinePar_t;
using CosmicSeedPar_t = Run1BAnaStructs::CosmicSeedPar_t;
using TimeClusterPar_t = Run1BAnaStructs::TimeClusterPar_t;
using SimPar_t = Run1BAnaStructs::SimPar_t;
using ChargedTrack_t = Run1BAnaStructs::ChargedTrack_t;
using Photon_t = Run1BAnaStructs::Photon_t;

namespace mu2e
{
  class Run1BAna : public art::EDAnalyzer
  {
  public:

    //--------------------------------------------------------------------------------------
    // Input config
    //--------------------------------------------------------------------------------------

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>      simCol           { Name("simCollection")          , Comment("simCollection")                       };
      fhicl::Atom<art::InputTag>      primary          { Name("primaryParticle")        , Comment("primaryParticle")                     };
      fhicl::Atom<art::InputTag>      digiCol          { Name("strawDigiMCCollection")  , Comment("strawDigiMCCollection")               };
      fhicl::Atom<art::InputTag>      comboCol         { Name("ComboHitCollection")     , Comment("ComboHitCollection")                  };
      fhicl::Atom<art::InputTag>      caloHitCol       { Name("CaloHitCollection")      , Comment("CaloHitCollection")                   };
      fhicl::Atom<art::InputTag>      caloClusterCol   { Name("caloClusterCollection")  , Comment("caloClusterCollection")               };
      fhicl::Atom<art::InputTag>      caloClusterMCCol { Name("caloClusterMCCollection"), Comment("caloClusterMCCollection")             };
      fhicl::Atom<art::InputTag>      caloShowerSimCol { Name("caloShowerSimCollection"), Comment("CaloShowerSimCollection")             };
      fhicl::Atom<art::InputTag>      lineCol          { Name("LineCollection")         , Comment("Kinematic line collection")           };
      fhicl::Atom<art::InputTag>      cosmicSeedCol    { Name("cosmicSeedCollection")   , Comment("Cosmic seed collection")              };
      fhicl::Atom<art::InputTag>      timeClusterCol   { Name("timeClusterCollection")  , Comment("Time cluster collection")             };
      fhicl::Atom<std::string>        trigProcess      { Name("triggerProcess")         , Comment("Process name for the trigger results")};
      fhicl::Atom<art::InputTag>      pbi              { Name("PBI")                    , Comment("ProtonBunchIntensity tag")            };
      fhicl::Atom<art::InputTag>      genCounter       { Name("genCounter")             , Comment("Generator counter tag")              , "genCounter"};
      fhicl::Atom<double>             maxGenEnergy     { Name("maxGenEnergy")           , Comment("Cut on the maximum primary energy")  , -1.};
      fhicl::Atom<int>                debugLevel       { Name("debugLevel")             , Comment("debugLevel")                         , 0 };
    };


    // Histograms
    struct Hist_t {
      EventHist_t   event;
      ClusterHist_t cluster;
      LineHist_t    line;
      CosmicSeedHist_t cosmic_seed;
      TimeClusterHist_t time_cluster;
      SimHist_t     sim;
    };

    //--------------------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------------------

    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit Run1BAna(const Parameters& conf);
    virtual void analyze(const art::Event& event) override;
    virtual void beginJob();
    virtual void endJob();
    virtual void beginRun(const art::Run&   run   );
    virtual void endRun(const art::Run& run ) override;

  private:

    void bookClusterHistograms    (ClusterHist_t*     Hist, const int index, const char* name);
    void bookEventHistograms      (EventHist_t*       Hist, const int index, const char* name);
    void bookLineHistograms       (LineHist_t*        Hist, const int index, const char* name);
    void bookCosmicSeedHistograms (CosmicSeedHist_t*  Hist, const int index, const char* name);
    void bookTimeClusterHistograms(TimeClusterHist_t* Hist, const int index, const char* name);
    void bookSimHistograms        (SimHist_t*         Hist, const int index, const char* name);
    void bookHistograms           (                         const int index, const char* name);

    void fillClusterHistograms     (ClusterHist_t*     Hist, double Weight = 1.);
    void fillLineHistograms        (LineHist_t*        Hist, double Weight = 1.);
    void fillCosmicSeedHistograms  (CosmicSeedHist_t*  Hist, double Weight = 1.);
    void fillTimeClusterHistograms (TimeClusterHist_t* Hist, double Weight = 1.);
    void fillEventHistograms       (EventHist_t*       Hist, double Weight = 1.);
    void fillSimHistograms         (SimHist_t*         Hist, double Weight = 1.);
    void fillHistograms            (Hist_t*            Hist);

    bool isGoodCluster(const CaloCluster* cluster);
    bool isGoodLine(const KalSeed* seed);
    bool isGoodCosmicSeed(const CosmicTrackSeed* seed);
    bool isGoodTimeCluster(const TimeCluster* tc);
    CLHEP::HepLorentzVector lineAtCluster(const CaloCluster* cl, const KalSeed* seed);
    CLHEP::HepLorentzVector lineSeedAtCluster(const CaloCluster* cl, const CosmicTrackSeed* seed);
    void initClusterPar(ClusterPar_t& par, const CaloCluster* cluster);
    void initLinePar(LinePar_t& par, const KalSeed* line);
    void initCosmicSeedPar(CosmicSeedPar_t& par, const CosmicTrackSeed* seed);
    void initTimeClusterPar(TimeClusterPar_t& par, const TimeCluster* tc);
    void matchLineToCluster(ClusterPar_t& par, const KalSeedCollection* lines);
    void matchSeedToCluster(ClusterPar_t& par, const CosmicTrackSeedCollection* seeds);
    void matchTimeClusterToCluster(ClusterPar_t& par, const TimeClusterCollection* time_clusters);
    float getTotalEnergyDepositedBySim(const CaloClusterMC* mc, const SimParticle* sim);
    float getAverageTimeDepositedBySim(const CaloClusterMC* mc, const SimParticle* sim);
    float getAverageTimeDeposited(const CaloClusterMC* mc);
    void getShowerSimEnergyAndAvgTime(const CaloShowerSimCollection* col, const SimParticle* sim, float& energy, float& avgTime) const;
    CLHEP::Hep3Vector getCrystalPosition(const int crystalID) const;
    CLHEP::Hep3Vector getSimParticleHitPosition(const CaloShowerSimCollection* col, const SimParticle* sim) const;

    // check if sim 1 or any of its heritage matches sim 2
    bool isRelated(const SimParticle* sim1, const SimParticle* sim2) const {
      if(!sim1 || !sim2) return false;
      if(&(*sim1) == &(*sim2)) return true; // if they match, return true
      if(sim1->hasParent() && isRelated(&(*sim1->parent()), sim2)) return true;
      return false;
    }

    //--------------------------------------------------------------------------------------
    // Data
    //--------------------------------------------------------------------------------------

    art::InputTag  sim_tag_;
    art::InputTag  primary_tag_;
    art::InputTag  mc_digi_tag_;
    art::InputTag  combo_hits_tag_;
    art::InputTag  calo_hits_tag_;
    art::InputTag  calo_cluster_mc_tag_;
    art::InputTag  calo_shower_sim_tag_;
    art::InputTag  clusters_tag_;
    art::InputTag  line_tag_;
    art::InputTag  cosmic_seed_tag_;
    art::InputTag  time_cluster_tag_;
    art::InputTag  trig_tag_;
    art::InputTag  pbi_tag_;
    double         max_gen_energy_;
    int            debug_level_;

    enum {kMaxHists = 100};
    Hist_t* hist_[kMaxHists];
    TH1* hist_norm_;
    const SimParticleCollection*           sim_col_ = nullptr;
    const PrimaryParticle*                 primary_ = nullptr;
    const StrawDigiMCCollection*           mc_digi_col_ = nullptr;
    const ComboHitCollection*              combo_hit_col_ = nullptr;
    const CaloHitCollection*               calo_hit_col_ = nullptr;
    const CaloClusterCollection*           cluster_col_ = nullptr;
    const CaloClusterMCCollection*         calo_cluster_mc_col_ = nullptr;
    const CaloClusterMCTruthAssn*          calo_cluster_mc_assn_ = nullptr;
    const CaloShowerSimCollection*         calo_shower_sim_col_ = nullptr;
    const KalSeedCollection*               line_col_ = nullptr;
    const CosmicTrackSeedCollection*       cosmic_seed_col_ = nullptr;
    const KalLineAssns*                    line_seed_assn_ = nullptr;
    const TimeClusterCollection*           time_cluster_col_ = nullptr;
    const TriggerResultsNavigator*         trig_nav_ = nullptr;
    const art::Event*                      event_ = nullptr;
    const Calorimeter*                     calorimeter_ = nullptr;
    EventPar_t                             evt_par_;
    ClusterPar_t                           cluster_par_;
    LinePar_t                              line_par_;
    CosmicSeedPar_t                        cosmic_seed_par_;
    TimeClusterPar_t                       time_cluster_par_;
    SimPar_t                               sim_par_;
    Photon_t                               photon_;
    std::map<unsigned, SimUtils::Sim_t>    sim_info_;
    mu2e::StopWatch*                       watch_; // for timing

    unsigned  long                         nevt_, ngen_;
  };

  //--------------------------------------------------------------------------------------
  Run1BAna::Run1BAna(const Parameters& config):
    art::EDAnalyzer{config}
    , sim_tag_            (config().simCol())
    , primary_tag_        (config().primary())
    , mc_digi_tag_        (config().digiCol())
    , combo_hits_tag_     (config().comboCol())
    , calo_hits_tag_      (config().caloHitCol())
    , calo_cluster_mc_tag_(config().caloClusterMCCol())
    , calo_shower_sim_tag_(config().caloShowerSimCol())
    , clusters_tag_       (config().caloClusterCol())
    , line_tag_           (config().lineCol())
    , cosmic_seed_tag_    (config().cosmicSeedCol())
    , time_cluster_tag_   (config().timeClusterCol())
    , trig_tag_           ("TriggerResults::" + config().trigProcess())
    , pbi_tag_            (config().pbi())
    , max_gen_energy_     (config().maxGenEnergy())
    , debug_level_        (config().debugLevel())
    , nevt_              (0)
  {
    // Normalization histogram
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir("data", "Data");
    hist_norm_ = dir.make<TH1D>("norm", "Normalization counts", 1, 0., 1.);

    // Analysis histograms
    for(int i = 0; i < kMaxHists; ++i) hist_[i] = nullptr;
    bookHistograms(0, "all_events");
    bookHistograms(1, "all_clusters");
    bookHistograms(2, "70MeV_clusters");
    bookHistograms(3, "70MeV_disk0");
    bookHistograms(4, "70MeV_disk1");
    bookHistograms(5, "cluster_id");
    bookHistograms(6, "time_30MeV_clusters");
    bookHistograms(7, "time_50MeV_clusters");
    bookHistograms(8, "10_sim_hits");
    bookHistograms(9, "10_sim_hits_70MeV");

    bookHistograms(10, "max_70MeV");
    bookHistograms(11, "max_70MeV_line");
    bookHistograms(12, "max_70MeV_noline");

    bookHistograms(20, "cluster_ID");
    bookHistograms(21, "cluster_ID_line");
    bookHistograms(22, "cluster_ID_noline");
    bookHistograms(23, "cluster_ID_80MeV");

    bookHistograms(25, "cluster_sig_ID");

    bookHistograms(30, "merged_clusters");
    bookHistograms(31, "merged_clusters_10MeV");

    bookHistograms(40, "unmerged_clusters");
    bookHistograms(41, "unmerged_clusters_30MeV");

    bookHistograms(60, "photon_");
    bookHistograms(61, "photon_no_weight");
    bookHistograms(62, "photon_no_tcl");
    bookHistograms(63, "photon_tcl_50ns");
    bookHistograms(64, "photon_tcl_50ns_100ns");
    bookHistograms(65, "photon_no_csm");
    bookHistograms(66, "photon_no_line");
    bookHistograms(67, "photon_no_tcl_r500");
    bookHistograms(68, "photon_no_tcl_r550");

    bookHistograms(80, "all_lines");
    bookHistograms(81, "lines_cluster");
    bookHistograms(82, "lines_cluster_70MeV");

    bookHistograms(90, "all_cosmic_seeds");
    bookHistograms(91, "cosmic_seed_ID");
    bookHistograms(92, "cosmic_seed_cls");
    bookHistograms(93, "cosmic_seed_cls_70MeV");

    bookHistograms(95, "all_time_clusters");
    bookHistograms(96, "time_cluster_most_sim_hits");


    // For timing info
    watch_ = new mu2e::StopWatch();
    watch_->Calibrate();
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::beginJob() {
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::beginRun(const art::Run & run) {
    // Get the calorimeter geometry
    const mu2e::GeomHandle<mu2e::Calorimeter> cal;
    calorimeter_ = &*cal;
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookEventHistograms(EventHist_t* Hist, const int index, const char* title) {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(std::format("evt_{}", index), title);
    Hist->npot           = dir.make<TH1F>("npot"          , "N(POT);"                     , 1000,    0.,   1.e8);
    Hist->n_mc_digis     = dir.make<TH1F>("nmc_digis"     , "N(MC digis);"                , 1000,    0., 20000.);
    Hist->ncombo_hits    = dir.make<TH1F>("ncombo_hits"   , "N(combo hits);"              ,  500,    0., 20000.);
    Hist->ncalo_hits     = dir.make<TH1F>("ncalo_hits"    , "N(calo hits);"               ,  500,    0.,  5000.);
    Hist->nclusters      = dir.make<TH1D>("nclusters"     , "N(calo clusters);"           ,  200,    0.,   200.);
    Hist->ngood_clusters = dir.make<TH1D>("ngood_clusters", "N(calo clusters | ID);"      ,  200,    0.,   200.);
    Hist->nlines         = dir.make<TH1D>("nlines"        , "N(Kinematic lines);"         ,  100,    0.,   100.);
    Hist->ngood_lines    = dir.make<TH1D>("ngood_lines"   , "N(Kinematic lines | ID);"    ,  100,    0.,   100.);
    Hist->ncosmic_seeds  = dir.make<TH1D>("ncosmic_seeds"   , "N(Cosmic seeds);"          ,  100,    0.,   100.);
    Hist->ngood_cosmic_seeds = dir.make<TH1D>("ngood_cosmic_seeds"   , "N(Cosmic seeds | ID);"    ,  100,    0.,   100.);
    Hist->ntime_clusters = dir.make<TH1D>("ntime_clusters", "N(time clusters);"           ,  100,    0.,   100.);
    Hist->ngood_time_clusters = dir.make<TH1D>("ngood_time_clusters", "N(time clusters | ID);"    ,  100,    0.,   100.);

    Hist->trig_bits     = dir.make<TH1D>("trig_bits"    , "Trigger bits;"               , 1000,    0.,   1000);
    Hist->trig_paths    = dir.make<TH1D>("trig_paths"   , "Trigger paths;"              ,  100,    0.,    100);
    Hist->sim_dr_dt = dir.make<TH2F>("sim_dr_dt", "Sim distance vs #Delta t;#Deltat (ns);#Delta r (mm)", 200, -100., 100., 100, 0., 200.);
    Hist->hit_x_y = dir.make<TH2F>("hit_x_y", "Sim hit positions in X-Y plane;X (mm);Y (mm)", 100, -800., 800., 100, -800., 800.);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookClusterHistograms(ClusterHist_t* Hist, const int index, const char* title) {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(std::format("cls_{}", index), title);
    Hist->energy        = dir.make<TH1F>("energy"  , "Cluster energy;Energy (MeV);"      ,  300,    0.,    300.);
    Hist->mc_edep       = dir.make<TH1F>("mc_edep" , "Cluster MC energy;MC Energy (MeV);",  300,    0.,    300.);
    Hist->time          = dir.make<TH1F>("time"    , "Cluster time;Time (ns);"           ,  100,    0.,   2000.);
    Hist->radius        = dir.make<TH1F>("radius"  , "Cluster radial position;R (mm);"   ,  100,    0.,    700.);
    Hist->ncr           = dir.make<TH1D>("ncr"     , "Number of crystals;N(Crystals);"   ,   20,    0.,     20.);
    Hist->disk          = dir.make<TH1D>("disk"    , "Disk ID;Disk ID;"                  ,    3,   -1.,      2.);
    Hist->energy_per_crystal     = dir.make<TH1F>("energy_per_crystal", "Cluster crystal energy;Energy (MeV);", 300, 0., 300.);
    Hist->frac_first_crystal     = dir.make<TH1F>("frac_first_crystal", "First crystal energy / cluster energy;Frac;", 101, 0., 1.01);
    Hist->frac_first_two_crystals= dir.make<TH1F>("frac_first_two_crystals", "First two crystals energy / cluster energy;Frac;", 101, 0., 1.01);
    Hist->second_moment = dir.make<TH1F>("second_moment", "Cluster second moment;Second moment;", 200, 0., 1.e6);
    Hist->e1            = dir.make<TH1F>("e1", "Cluster E1;E1 (MeV);", 300, 0., 300.);
    Hist->e2            = dir.make<TH1F>("e2", "Cluster E2;E2 (MeV);", 300, 0., 300.);
    Hist->e9            = dir.make<TH1F>("e9", "Cluster E9;E9 (MeV);", 300, 0., 300.);
    Hist->e25           = dir.make<TH1F>("e25", "Cluster E25;E25 (MeV);", 300, 0., 300.);
    Hist->t_var         = dir.make<TH1F>("t_var", "Cluster time variance;#sigma_{t}^{2} (ns^{2});", 200, 0., 10.);
    Hist->line_dt       = dir.make<TH1F>("line_dt" , "Line-cluster time diff;dt (ns);"   ,  100, -200.,    200.);
    Hist->line_dr       = dir.make<TH1F>("line_dr" , "Line-cluster distance;dr (mm);"    ,  100,    0.,    500.);
    Hist->nmatched_lines= dir.make<TH1D>("nmatched_lines", "N(lines) matched to the cluster;N(lines)",  40, 0.,    40.);
    Hist->nfit_matched_lines= dir.make<TH1D>("nfit_matched_lines", "N(lines) fit matched to the cluster;N(lines)",  40, 0.,    40.);

    Hist->time_cluster_dt = dir.make<TH1F>("time_cluster_dt", "Time cluster-cluster time diff;dt (ns);",  100, -200.,    200.);
    Hist->time_cluster_dr = dir.make<TH1F>("time_cluster_dr", "Time cluster-cluster distance;dr (mm);",  100,    0.,    500.);
    Hist->nmatched_time_clusters= dir.make<TH1D>("nmatched_time_clusters", "N(time clusters) matched to the cluster;N(time clusters)",  40, 0.,    40.);
    Hist->nfit_matched_time_clusters= dir.make<TH1D>("nfit_matched_time_clusters", "N(time clusters) fit matched to the cluster;N(time clusters)",  40, 0.,    40.);

    Hist->esum          = dir.make<TH1F>("esum_cls"      , "Energy sum between two clusters;E_{sum} (MeV);"    ,  300,    0.,    300.);
    Hist->dt            = dir.make<TH1F>("dt_cls"        , "Time difference between two clusters;#Deltat (ns);",  300,    0.,    300.);
    Hist->dr            = dir.make<TH1F>("dr_cls"        , "Position difference between two clusters;#Deltar (mm);",  300, 0.,   300.);
    Hist->dr_vs_dt      = dir.make<TH2F>("dr_vs_dt_cls"  , "Position and time difference between two clusters;#Deltat (ns);#Deltar (mm)",  100, -100., 100., 100, 0., 200.);

    // MC info
    Hist->pdg           = dir.make<TH1D>("pdg"       , "Primary sim PDG ID;PDG ID;"                          , 2500, -200.,  2300.);
    Hist->energy_sim    = dir.make<TH1F>("energy_sim", "Primary sim total energy deposited;Energy (MeV);"    ,  300,    0.,    300.);
    Hist->energy_ratio  = dir.make<TH1F>("energy_ratio", "Energy ratio (primary/total);E_sim/E_total;"       ,  110,    0.,    1.1);
    Hist->pdg2          = dir.make<TH1D>("pdg2"      , "Secondary sim PDG ID;PDG ID;"                        , 2500, -200.,  2300.);
    Hist->energy_sim2   = dir.make<TH1F>("energy_sim2", "Secondary sim total energy deposited;Energy (MeV);"  ,  300,    0.,    300.);
    Hist->energy_ratio2 = dir.make<TH1F>("energy_ratio2", "Energy ratio (secondary/total);E_sim/E_total;"     ,  110,    0.,    1.1);
    Hist->sim_dt        = dir.make<TH1F>("sim_dt", "Primary-secondary sim time diff;#Delta t (ns);", 200, -50., 50.);
    Hist->MC_energy_diff= dir.make<TH1F>("energy_res", "Reconstructed - MC cluster energy;#Delta E (MeV);", 300, -15., 15.);
    Hist->MC_time_diff  = dir.make<TH1F>("time_res", "Reconstructed - MC cluster time;#Delta t (ns);", 200, -20., 20.);
    Hist->energy_vs_gen_energy = dir.make<TH2F>("energy_vs_gen_energy", "Cluster energy vs. generated energy;E_{gen} (MeV);E_{cluster} (MeV)",  60, 50., 110., 80, 40., 120.);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookLineHistograms(LineHist_t* Hist, const int index, const char* title) {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(std::format("line_{}", index), title);
    Hist->chi2       = dir.make<TH1F>("chi2",       "Line fit chi2;chi2;"                  , 100,   0.,     10.);
    Hist->nhits      = dir.make<TH1F>("nhits",      "Line nhits;N(straw hits);"            , 100,   0.,    100.);
    Hist->nplanes    = dir.make<TH1F>("nplanes",    "Line nplanes;N(planes);"              ,  50,   0.,     50.);
    Hist->nstereo    = dir.make<TH1F>("nstereo",    "Line nstereo;N(stereo panels);"       ,  15,   0.,     15.);
    Hist->d0         = dir.make<TH1F>("d0",         "Line d0;d0 (mm);"                     , 200, -400.,    400.);
    Hist->tdip       = dir.make<TH1F>("tdip",       "Line tan(dip);tan(dip);"              , 100, -10.,     10.);
    Hist->cos        = dir.make<TH1F>("cos",        "Line cos(theta);cos(theta);"          , 100,   0.,      1.);
    Hist->z0         = dir.make<TH1F>("z0",         "Line z0;z0 (mm);"                     , 100,-5000.,   5000.);
    Hist->t0         = dir.make<TH1F>("t0",         "Line t0;t0 (ns);"                     , 200,    0.,   2000.);
    Hist->phi0       = dir.make<TH1F>("phi0",       "Line phi0;phi0 (rad);"                , 100,  -3.15,   3.15);
    Hist->cl_energy  = dir.make<TH1F>("cl_energy",  "Cluster energy from line match;E (MeV);", 300, 0., 300.);
    Hist->cl_dt      = dir.make<TH1F>("cl_dt",      "Line-cluster time diff;dt (ns);"     , 200, -200.,    200.);
    Hist->cl_dr      = dir.make<TH1F>("cl_dr",      "Line-cluster distance;dr (mm);"      , 150,   0.,    500.);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookCosmicSeedHistograms(CosmicSeedHist_t* Hist, const int index, const char* title) {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(std::format("csms_{}", index), title);
    Hist->chi2       = dir.make<TH1F>("chi2",       "Cosmic seed fit chi2;chi2;"                  , 100,   0.,     10.);
    Hist->nhits      = dir.make<TH1F>("nhits",      "Cosmic seed nhits;N(straw hits);"            , 100,   0.,    100.);
    Hist->d0         = dir.make<TH1F>("d0",         "Cosmic seed d0;d0 (mm);"                     , 200, -400.,    400.);
    Hist->tdip       = dir.make<TH1F>("tdip",       "Cosmic seed tan(dip);tan(dip);"              , 100, -10.,     10.);
    Hist->cos        = dir.make<TH1F>("cos",        "Cosmic seed cos(theta);cos(theta);"          , 100,   0.,      1.);
    Hist->z0         = dir.make<TH1F>("z0",         "Cosmic seed z0;z0 (mm);"                     , 100,-5000.,   5000.);
    Hist->t0         = dir.make<TH1F>("t0",         "Cosmic seed t0;t0 (ns);"                     , 200,    0.,   2000.);
    Hist->phi0       = dir.make<TH1F>("phi0",       "Cosmic seed phi0;phi0 (rad);"                , 100,  -3.15,   3.15);
    Hist->A0         = dir.make<TH1F>("A0",         "Cosmic seed A0;A0 (mm);"                     , 200, -4000.,  4000.);
    Hist->B0         = dir.make<TH1F>("B0",         "Cosmic seed B0;B0 (mm);"                     , 200, -4000.,  4000.);
    Hist->A1         = dir.make<TH1F>("A1",         "Cosmic seed A1;A1 (mm);"                     , 200, -4000.,  4000.);
    Hist->B1         = dir.make<TH1F>("B1",         "Cosmic seed B1;B1 (mm);"                     , 200, -4000.,  4000.);
    Hist->cl_energy  = dir.make<TH1F>("cl_energy",  "Cluster energy from cosmic seed match;E (MeV);", 300, 0., 300.);
    Hist->cl_time    = dir.make<TH1F>("cl_time",    "Cluster time from cosmic seed match;time (ns);", 200, 0., 2000.);
    Hist->cl_disk    = dir.make<TH1D>("cl_disk",    "Cluster disk from cosmic seed match;disk ID;"  , 5, -1., 4.);
    Hist->cl_dt      = dir.make<TH1F>("cl_dt",      "Cosmic seed-cluster time diff;dt (ns);"     , 200, -200.,    200.);
    Hist->cl_dr      = dir.make<TH1F>("cl_dr",      "Cosmic seed-cluster distance;dr (mm);"      , 150,   0.,    500.);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookTimeClusterHistograms(TimeClusterHist_t* Hist, const int index, const char* title) {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(std::format("tcls_{}", index), title);
    Hist->nhits          = dir.make<TH1F>("nhits",          "Time cluster nhits;N(hits);"                    , 200,    0.,    200.);
    Hist->nstraw_hits    = dir.make<TH1F>("nstraw_hits",    "Time cluster straw hits;N(straw hits);"         , 200,    0.,    200.);
    Hist->t0             = dir.make<TH1F>("t0",             "Time cluster t0;t0 (ns);"                       , 200,    0.,   2000.);
    Hist->t0err          = dir.make<TH1F>("t0err",          "Time cluster t0 error;t0 error (ns);"           , 100,    0.,     10.);
    Hist->z0             = dir.make<TH1F>("z0",             "Time cluster z0;z0 (mm);"                       , 100,-5000.,   5000.);
    Hist->phi0           = dir.make<TH1F>("phi0",           "Time cluster phi0;phi0;"                        , 100, -3.15,    3.15);
    Hist->cl_energy      = dir.make<TH1F>("cl_energy",      "Cluster energy from time cluster match;E (MeV);", 300,    0.,    300.);
    Hist->cl_time        = dir.make<TH1F>("cl_time",        "Cluster time from time cluster match;time (ns);", 200,    0.,   2000.);
    Hist->cl_disk        = dir.make<TH1D>("cl_disk",        "Cluster disk from time cluster match;disk ID;"  ,   5,   -1.,      4.);
    Hist->cl_dt          = dir.make<TH1F>("cl_dt",          "Time cluster-cluster time diff;dt (ns);"        , 200, -200.,    200.);
    Hist->cl_dr          = dir.make<TH1F>("cl_dr",          "Time cluster-cluster distance;dr (mm);"         , 150,    0.,    500.);
    Hist->n_primary_hits = dir.make<TH1F>("n_primary_hits", "Time cluster primary hits;N(primary hits);"     , 200,    0.,    200.);
    Hist->n_other_hits   = dir.make<TH1F>("n_other_hits",   "Time cluster other hits;N(other hits);"         , 200,    0.,    200.);
    Hist->purity         = dir.make<TH1F>("purity",         "Time cluster purity;Purity;"                    , 110,    0.,     1.1);
    Hist->efficiency     = dir.make<TH1F>("efficiency",     "Time cluster efficiency;Efficiency;"            , 110,    0.,     1.1);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookSimHistograms(SimHist_t* Hist, const int index, const char* title) {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(std::format("sim_{}", index), title);
    Hist->nhits        = dir.make<TH1D>("nhits"       , "Sim N(hits);N(hits);"          ,  100,    0.,   200.);
    Hist->pdg          = dir.make<TH1D>("pdg"         , "Sim PDG ID;ID;"                , 2500, -200.,  2300.);
    Hist->type         = dir.make<TH1D>("type"         , "Sim type;;"                    ,   10,    0.,    10.);
    Hist->parent_type  = dir.make<TH1D>("parent_type" , "Sim parent type;;"             ,   10,    0.,    10.);
    Hist->origin_type  = dir.make<TH1D>("origin_type" , "Sim origin type;;"             ,   10,    0.,    10.);
    Hist->energy_start = dir.make<TH1D>("energy_start", "Sim start energy;;"            ,  150,    0.,   150.);
    Hist->start_x_y    = dir.make<TH2F>("start_x_y"   , "Sim origin;x (mm);y (mm)"      ,   80, -200., 200., 80, -200., 200.);
    Hist->energy_vs_trig_path = dir.make<TH2F>("energy_vs_trig_path", "Trigger path vs. gen energy;;E_{gen} (MeV)", 20., 0., 20., 60, 50., 110.);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookHistograms(const int index, const char* title) {
    if(index >= kMaxHists) throw std::runtime_error("Too many histograms!");
    hist_[index] = new Hist_t;
    auto Hist = hist_[index];
    bookEventHistograms      (&Hist->event       , index, title);
    bookClusterHistograms    (&Hist->cluster     , index, title);
    bookLineHistograms       (&Hist->line        , index, title);
    bookCosmicSeedHistograms (&Hist->cosmic_seed , index, title);
    bookTimeClusterHistograms(&Hist->time_cluster, index, title);
    bookSimHistograms        (&Hist->sim         , index, title);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillClusterHistograms(ClusterHist_t* Hist, double Weight) {
    if(!Hist) return;
    auto Cluster = cluster_par_.cluster;
    auto Line = cluster_par_.line;
    auto TimeCluster = cluster_par_.time_cluster;
    if(!Cluster) return;

    Hist->energy->Fill(Cluster->energyDep(), Weight);
    if(cluster_par_.mc) Hist->mc_edep->Fill(cluster_par_.mc->totalEnergyDep(), Weight);
    Hist->time->Fill(Cluster->time(), Weight);
    Hist->radius->Fill(cluster_par_.r, Weight);
    Hist->ncr->Fill(Cluster->caloHitsPtrVector().size(), Weight);
    Hist->disk->Fill(Cluster->diskID(), Weight);
    Hist->second_moment->Fill(cluster_par_.second_moment, Weight);
    Hist->e1->Fill(cluster_par_.e1, Weight);
    Hist->e2->Fill(cluster_par_.e2, Weight);
    Hist->e9->Fill(cluster_par_.e9, Weight);
    Hist->e25->Fill(cluster_par_.e25, Weight);
    Hist->t_var->Fill(cluster_par_.t_var, Weight);
    if(Line) {
      const auto line_pos = lineAtCluster(Cluster, Line);
      const auto cl_pos = Cluster->cog3Vector();
      const double dx = line_pos.x() - cl_pos.x();
      const double dy = line_pos.y() - cl_pos.y();
      const double dr = std::sqrt(dx*dx + dy*dy);
      const double dt = line_pos.t() - Cluster->time();

      Hist->line_dt->Fill(dt, Weight);
      Hist->line_dr->Fill(dr, Weight);
    }
    Hist->nmatched_lines->Fill(cluster_par_.nmatched_lines, Weight);
    Hist->nfit_matched_lines->Fill(cluster_par_.nfit_matched_lines, Weight);

    if(TimeCluster) {
      const auto tc_pos = TimeCluster->position();
      const auto cl_pos = Cluster->cog3Vector();
      const double dx = tc_pos.x() - cl_pos.x();
      const double dy = tc_pos.y() - cl_pos.y();
      const double dr = std::sqrt(dx*dx + dy*dy);
      const double dt = TimeCluster->t0().t0() - Cluster->time();

      Hist->time_cluster_dt->Fill(dt, Weight);
      Hist->time_cluster_dr->Fill(dr, Weight);
    }
    Hist->nmatched_time_clusters->Fill(cluster_par_.nmatched_time_clusters, Weight);
    Hist->nfit_matched_time_clusters->Fill(cluster_par_.nfit_matched_time_clusters, Weight);

    // Per-cluster per-crystal and fraction metrics
    const auto& chptrs = Cluster->caloHitsPtrVector();
    const size_t ncr = chptrs.size();
    const float per_crystal = Cluster->energyDep() / static_cast<float>(ncr);
    Hist->energy_per_crystal->Fill(per_crystal, Weight);
    Hist->frac_first_crystal->Fill(cluster_par_.frac_1(), Weight);
    Hist->frac_first_two_crystals->Fill(cluster_par_.frac_2(), Weight);

    for(const auto& cls : *cluster_col_) {
      if(&cls == &(*Cluster)) continue;
      if(cls.energyDep() < 30.) continue; // only look at other clusters with >30 MeV energy to avoid too much background
      if(cls.diskID() != Cluster->diskID()) continue; // only look at clusters on the same disk
      const float dx = (Cluster->cog3Vector() - cls.cog3Vector()).x();
      const float dy = (Cluster->cog3Vector() - cls.cog3Vector()).y();
      const float dr = std::sqrt(dx*dx + dy*dy);
      const float dt = Cluster->time() - cls.time();
      Hist->esum      ->Fill(Cluster->energyDep() + cls.energyDep(), Weight);
      Hist->dt        ->Fill(dt, Weight);
      Hist->dr        ->Fill(dr, Weight);
      Hist->dr_vs_dt  ->Fill(dt, dr, Weight);
    }

    // MC information
    if(cluster_par_.mc) {
      const float mc_energy = cluster_par_.mc->totalEnergyDep();
      Hist->MC_energy_diff->Fill(Cluster->energyDep() - mc_energy, Weight);
      const float mc_time = getAverageTimeDeposited(cluster_par_.mc);
      Hist->MC_time_diff->Fill(Cluster->time() - mc_time, Weight);
      const auto sim_1 = cluster_par_.primary_sim;
      const auto sim_2 = cluster_par_.secondary_sim;
      const float sim_edep_1 = (sim_1) ? getTotalEnergyDepositedBySim(cluster_par_.mc, sim_1) : 0.;
      const float sim_edep_2 = (sim_2) ? getTotalEnergyDepositedBySim(cluster_par_.mc, sim_2) : 0.;
      const float sim_time_1 = (sim_1) ? getAverageTimeDepositedBySim(cluster_par_.mc, sim_1) : 0.;
      const float sim_time_2 = (sim_2) ? getAverageTimeDepositedBySim(cluster_par_.mc, sim_2) : 0.;
      const float sim_pdg_1  = (sim_1) ? int(sim_1->pdgId()) : 0.;
      const float sim_pdg_2  = (sim_2) ? int(sim_2->pdgId()) : 0.;

      Hist->pdg->Fill(sim_pdg_1, Weight);
      Hist->energy_sim->Fill(sim_edep_1, Weight);
      Hist->energy_ratio->Fill(sim_edep_1 / mc_energy, Weight);

      Hist->pdg2->Fill(sim_pdg_2, Weight);
      Hist->energy_sim2->Fill(sim_edep_2, Weight);
      Hist->energy_ratio2->Fill(sim_edep_2 / mc_energy, Weight);

      if((sim_edep_1 + sim_edep_2)/mc_energy > 1.01) {
        std::cout << "[Run1BAna::" << __func__ << "] " << event_->id()
                  << ": Unphysical energy fractions!" << std::endl
                  << " Cluster: E = " << Cluster->energyDep() << " Time = " << Cluster->time() << std::endl
                  << "  Sim 1: PDG = " << sim_pdg_1 << " E = " << sim_edep_1 << " Time = " << sim_time_1
                  << " Frac = " << sim_edep_1 / mc_energy << std::endl
                  << "  Sim 2: PDG = " << sim_pdg_2 << " E = " << sim_edep_2 << " Time = " << sim_time_2
                  << " Frac = " << sim_edep_2 / mc_energy << std::endl;
      }

      if(sim_1 && sim_2) { // only fill if both are found
        Hist->sim_dt->Fill(sim_time_1 - sim_time_2, Weight);
      }
    }

    const float e_gen = (cluster_par_.primary_sim) ? cluster_par_.primary_sim->startMomentum().e() : 0.;
    Hist->energy_vs_gen_energy->Fill(e_gen, Cluster->energyDep(), Weight);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillLineHistograms(LineHist_t* Hist, double Weight) {
    if(!Hist) return;
    auto Line = line_par_.line;
    if(!Line) return;
    auto Cluster = (Line->hasCaloCluster()) ? &(*Line->caloCluster()) : nullptr;

    // Extract line parameters from KalSeed
    double t0;
    try {
      Line->t0Segment(t0);
    } catch(...) { return; }

    auto t0seg = Line->t0Segment(t0);
    if(t0seg == Line->segments().end()) return;
    auto momvec = t0seg->momentum3();
    auto posvec = t0seg->position3();
    if(std::abs(posvec.x()) > 1.e10) return; // bad line
    double theta = momvec.Theta();
    double phi = momvec.Phi();
    double td = 1.0 / tan(theta);

    double d0;
    try {
      auto kltraj = t0seg->kinematicLine();
      d0 = kltraj.d0();
    } catch (...) { return; }

    // Count planes and stereo panels
    std::set<unsigned> stcount;
    std::set<unsigned> pcount;
    for(const auto& hit : Line->hits()) {
      if(hit._flag.hasAllProperties(StrawHitFlag::active)) {
        stcount.insert(hit._sid.stereoPanel());
        pcount.insert(hit._sid.plane());
      }
    }

    unsigned nactive = Line->nHits(true);

    // Fill histograms
    Hist->chi2->Fill(Line->chisquared() / Line->nDOF(), Weight);
    Hist->nhits->Fill(nactive, Weight);
    Hist->nplanes->Fill(pcount.size(), Weight);
    Hist->nstereo->Fill(stcount.size(), Weight);
    Hist->d0->Fill(d0, Weight);
    Hist->tdip->Fill(td, Weight);
    Hist->cos->Fill(std::cos(theta), Weight);
    Hist->z0->Fill(posvec.Z(), Weight);
    Hist->t0->Fill(t0, Weight);
    Hist->phi0->Fill(phi, Weight);

    // Fill cluster-matched histograms if available
    if(Cluster && Cluster->energyDep() > 1.) {
      const auto line_pos = lineAtCluster(Cluster, Line);
      const auto cl_pos = Cluster->cog3Vector();
      const double dx = line_pos.x() - cl_pos.x();
      const double dy = line_pos.y() - cl_pos.y();
      const double dr = std::sqrt(dx*dx + dy*dy);
      const double dt = line_pos.t() - Cluster->time();


      Hist->cl_energy->Fill(Cluster->energyDep(), Weight);
      Hist->cl_dt->Fill(dt, Weight);
      Hist->cl_dr->Fill(dr, Weight);
    }
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillCosmicSeedHistograms(CosmicSeedHist_t* Hist, double Weight) {
    if(!Hist) return;
    auto CosmicSeed = cosmic_seed_par_.seed;
    if(!CosmicSeed) return;
    const auto& Track = CosmicSeed->track();
    auto Cluster = (CosmicSeed->hasCaloCluster()) ? &(*CosmicSeed->caloCluster()) : nullptr;
    Hist->nhits->Fill(CosmicSeed->hits().size(), Weight);
    Hist->t0->Fill(CosmicSeed->t0().t0(), Weight);
    Hist->d0->Fill(Track.d0(), Weight);
    Hist->phi0->Fill(Track.phi0(), Weight);
    Hist->z0->Fill(Track.z0(), Weight);
    Hist->cos->Fill(Track.cost(), Weight);
    Hist->A0->Fill(Track.FitParams.A0, Weight);
    Hist->B0->Fill(Track.FitParams.B0, Weight);
    Hist->A1->Fill(Track.FitParams.A1, Weight);
    Hist->B1->Fill(Track.FitParams.B1, Weight);
    // Hist->t0err->Fill(CosmicSeed->t0().t0Err(), Weight);
    // Hist->z0->Fill(CosmicSeed->position().z(), Weight);
    // Hist->phi0->Fill(CosmicSeed->position().phi(), Weight);
    if(Cluster) {
      Hist->cl_energy->Fill(Cluster->energyDep(), Weight);
      Hist->cl_time->Fill(Cluster->time(), Weight);
      Hist->cl_disk->Fill(Cluster->diskID(), Weight);
      const double dt = CosmicSeed->t0().t0() - Cluster->time();
      // Hist->cl_dr     ->Fill(dr, Weight);
      Hist->cl_dt     ->Fill(dt, Weight);
    }
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillTimeClusterHistograms(TimeClusterHist_t* Hist, double Weight) {
    if(!Hist) return;
    auto TimeCluster = time_cluster_par_.tc;
    if(!TimeCluster) return;
    auto Cluster = (TimeCluster->hasCaloCluster()) ? &(*TimeCluster->caloCluster()) : nullptr;

    Hist->nhits->Fill(TimeCluster->hits().size(), Weight);
    Hist->nstraw_hits->Fill(TimeCluster->nStrawHits(), Weight);
    Hist->t0->Fill(TimeCluster->t0().t0(), Weight);
    Hist->t0err->Fill(TimeCluster->t0().t0Err(), Weight);
    Hist->z0->Fill(TimeCluster->position().z(), Weight);
    Hist->phi0->Fill(TimeCluster->position().phi(), Weight);
    if(Cluster) {
      Hist->cl_energy->Fill(Cluster->energyDep(), Weight);
      Hist->cl_time->Fill(Cluster->time(), Weight);
      Hist->cl_disk->Fill(Cluster->diskID(), Weight);
      const auto tc_pos = TimeCluster->position();
      const auto cl_pos = Cluster->cog3Vector();
      const double dx = tc_pos.x() - cl_pos.x();
      const double dy = tc_pos.y() - cl_pos.y();
      const double dr = std::sqrt(dx*dx + dy*dy);
      const double dt = TimeCluster->t0().t0() - Cluster->time();
      Hist->cl_dr->Fill(dr, Weight);
      Hist->cl_dt->Fill(dt, Weight);
    }
    Hist->n_primary_hits->Fill(time_cluster_par_.n_primary_hits, Weight);
    Hist->n_other_hits->Fill(time_cluster_par_.n_other_hits, Weight);
    Hist->purity->Fill(time_cluster_par_.purity(), Weight);
    Hist->efficiency->Fill(time_cluster_par_.efficiency(), Weight);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillEventHistograms(EventHist_t* Hist, double Weight) {
    if(!Hist) return;
    Hist->npot           ->Fill(evt_par_.npot, Weight);
    Hist->n_mc_digis     ->Fill((mc_digi_col_) ? mc_digi_col_->size() : 0, Weight);
    Hist->ncombo_hits    ->Fill((combo_hit_col_) ? combo_hit_col_->size() : 0, Weight);
    Hist->ncalo_hits     ->Fill((calo_hit_col_) ? calo_hit_col_->size() : 0, Weight);
    Hist->nclusters      ->Fill((cluster_col_) ? cluster_col_->size() : 0, Weight);
    Hist->ngood_clusters ->Fill(evt_par_.n_good_clusters, Weight);
    Hist->nlines         ->Fill((line_col_) ? line_col_->size() : 0, Weight);
    Hist->ngood_lines    ->Fill(evt_par_.n_good_lines, Weight);
    Hist->ncosmic_seeds  ->Fill((cosmic_seed_col_) ? cosmic_seed_col_->size() : 0, Weight);
    Hist->ngood_cosmic_seeds->Fill(evt_par_.n_good_cosmic_seeds, Weight);
    Hist->ntime_clusters ->Fill((time_cluster_col_) ? time_cluster_col_->size() : 0, Weight);
    Hist->ngood_time_clusters ->Fill(evt_par_.n_good_time_clusters, Weight);

    // Trigger information
    for (size_t index = 0; index < trig_nav_->getTrigPaths().size(); ++index) {
      const std::string path = trig_nav_->getTrigPathNameByIndex(index);
      if(trig_nav_->accepted(path)) {
        Hist->trig_bits ->Fill(trig_nav_->getTrigBitByName(path), Weight);
        Hist->trig_paths->Fill(path.c_str(), Weight);
      }
    }

    // Event-level primary sim distance/time difference 2D histogram
    if(primary_ && calo_shower_sim_col_ && primary_->primarySimParticles().size() > 1) {
      const auto& pvec = primary_->primarySimParticles();
      const SimParticle* psim0 = &(*pvec.front());
      const SimParticle* psim1 = &(*pvec.at(1));
      float e0 = 0.f, t0 = 0.f;
      float e1 = 0.f, t1 = 0.f;
      getShowerSimEnergyAndAvgTime(calo_shower_sim_col_, psim0, e0, t0);
      getShowerSimEnergyAndAvgTime(calo_shower_sim_col_, psim1, e1, t1);
      const auto pos0 = getSimParticleHitPosition(calo_shower_sim_col_, psim0);
      const auto pos1 = getSimParticleHitPosition(calo_shower_sim_col_, psim1);
      const double dx = pos0.x() - pos1.x();
      const double dy = pos0.y() - pos1.y();
      const double dz = pos0.z() - pos1.z();
      const double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
      if(e0 > 0.f && e1 > 0.f) {
        Hist->sim_dr_dt->Fill(t0 - t1, dist, Weight);
        Hist->hit_x_y->Fill(pos0.x() + 3904., pos0.y(), Weight);
        Hist->hit_x_y->Fill(pos1.x() + 3904., pos1.y(), Weight);
      }
    }
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillSimHistograms(SimHist_t* Hist, double Weight) {
    if(!Hist) return;
    const SimParticle* sim = sim_par_.sim;
    if(!sim) return;
    Hist->pdg         ->Fill(sim->pdgId(), Weight);
    Hist->type        ->Fill(SimUtils::getSimType(sim), Weight);
    Hist->parent_type ->Fill(SimUtils::getSimType((sim->hasParent()) ? &(*sim->parent()) : nullptr), Weight);
    Hist->origin_type ->Fill(SimUtils::getSimOriginType(sim), Weight);
    Hist->energy_start->Fill(sim->startMomentum().e(), Weight);
    Hist->nhits       ->Fill(sim_par_.nhits, Weight);
    Hist->start_x_y   ->Fill(sim->startPosition().x() + 3904., sim->startPosition().y(), Weight);

    // Trigger information
    for (size_t index = 0; index < trig_nav_->getTrigPaths().size(); ++index) {
      const std::string path = trig_nav_->getTrigPathNameByIndex(index);
      if(trig_nav_->accepted(path)) {
        Hist->energy_vs_trig_path->Fill(path.c_str(), sim->startMomentum().e(), Weight);
      }
    }
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillHistograms(Hist_t* Hist) {
    if(!Hist) return;
    const double weight = evt_par_.weight;
    watch_->SetTime("FillHistogram");
    fillEventHistograms(&Hist->event, weight);
    fillClusterHistograms(&Hist->cluster, weight);
    fillLineHistograms(&Hist->line, weight);
    fillCosmicSeedHistograms(&Hist->cosmic_seed, weight);
    fillTimeClusterHistograms(&Hist->time_cluster, weight);
    fillSimHistograms(&Hist->sim, weight);
    watch_->StopTime("FillHistogram");
  }

  //--------------------------------------------------------------------------------------
  bool Run1BAna::isGoodTimeCluster(const TimeCluster* tc) {
    if(!tc) return false;
    if(tc->nhits() < 3) return false; // must have at least 3 hits
    // if(tc->t0().t0() < 600. || tc->t0().t0() > 1650.) return false; // time window selection
    return true;
  }

  //--------------------------------------------------------------------------------------
  bool Run1BAna::isGoodCosmicSeed(const CosmicTrackSeed* seed) {
    if(!seed) return false;
    if(seed->hits().size() < 5) return false; // must have at least 5 hits
    // if(seed->t0().t0() < 550. || seed->t0().t0() > 1650.) return false; // time window selection
    return true; // passes all selections
  }

  //--------------------------------------------------------------------------------------
  bool Run1BAna::isGoodLine(const KalSeed* seed) {
    if(!seed) return false;
    if(seed->intersections().empty()) return false; // intersections list must be populated

    constexpr int min_nactive(4), min_nstereo(1), min_npanels(1);
    // constexpr int min_nactive(5), min_nstereo(2), min_npanels(3);
    constexpr double max_chi2_dof(10.);

    if(seed->nHits(true) < min_nactive) return false; // active hits
    if(seed->chisquared() / seed->nDOF() > max_chi2_dof) return false; // chi^2/dof
    // compute number of planes and unique stereo orientations
    std::set<unsigned> stcount;
    std::set<unsigned> pcount;
    for(auto const& hit : seed->hits()){
      if(hit._flag.hasAllProperties(StrawHitFlag::active)){
        stcount.insert(hit._sid.stereoPanel());
        pcount.insert(hit._sid.plane());
      }
    }
    if(stcount.size() < min_nstereo) return false;
    if(pcount .size() < min_npanels) return false;

    return true; // passes all selections
  }

  //--------------------------------------------------------------------------------------
  bool Run1BAna::isGoodCluster(const CaloCluster* cluster) {
    if(!cluster) return false;
    const int disk_id = cluster->diskID();
    if(disk_id != 0 && disk_id != 1) {
      std::cout << "[Run1BAna::" << __func__ << "] " << event_->id()
                << ": Bad cluster! Disk ID = " << disk_id << std::endl;
      return false;
    }

    return true; // passes all selections
  }

  //--------------------------------------------------------------------------------------
  CLHEP::HepLorentzVector Run1BAna::lineAtCluster(const CaloCluster* cl, const KalSeed* seed) {
    mu2e::GeomHandle<mu2e::Calorimeter> cal;
    return mu2e::Run1BAnaUtils::lineAtCluster(cl, seed, cal,debug_level_-1);
  }

  //--------------------------------------------------------------------------------------
  CLHEP::HepLorentzVector Run1BAna::lineSeedAtCluster(const CaloCluster* cl, const CosmicTrackSeed* seed) {
    mu2e::GeomHandle<mu2e::Calorimeter> cal;
    return Run1BAnaUtils::lineSeedAtCluster(cl, seed, cal, debug_level_-1);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::initClusterPar(ClusterPar_t& par, const CaloCluster* cluster) {
    par.init(cluster);
    matchLineToCluster(par, line_col_);
    matchSeedToCluster(par, cosmic_seed_col_);
    matchTimeClusterToCluster(par, time_cluster_col_);

    // Match reconstructed cluster to MC cluster using the CaloClusterMCTruthAssn
    if(!cluster) return;

    if(calo_cluster_mc_assn_) {
      for(const auto& ent : *calo_cluster_mc_assn_) {
        const art::Ptr<CaloCluster>& recoPtr = ent.first;
        const art::Ptr<CaloClusterMC>& mcPtr = ent.second;
        if(recoPtr.isNonnull() && mcPtr.isNonnull()) {
          if(&(*cluster) == &(*recoPtr)) {
            par.mc = &(*mcPtr);
            break;
          }
        }
      }
      if(!par.mc) {
        std::cerr << "Warning: no MC match found for cluster with energy " << par.cluster->energyDep() << " MeV at time " << par.cluster->time() << " ns" << std::endl;
      }
    } else if(calo_cluster_mc_col_) { // no associations, so just try by index alignment
      if(calo_cluster_mc_col_->size() != cluster_col_->size()) {
        std::cerr << "Warning: Cluster and MC Cluster collections not aligned: N(clusters) = " << cluster_col_->size()
                  << " N(MC clusters) = " << calo_cluster_mc_col_->size() << std::endl;
      } else {
        // find this cluster's index
        size_t index = 0;
        for(const auto& cl : *cluster_col_) {
          if(&cl == &(*cluster)) break;
          ++index;
        }
        if(index >= calo_cluster_mc_col_->size()) {
          std::cerr << "Warning: no MC cluster by index found for cluster with energy " << par.cluster->energyDep() << " MeV at time " << par.cluster->time() << " ns" << std::endl;
        } else {
          par.mc = &calo_cluster_mc_col_->at(index);
        }
      }
    }

    // Initialize MC info if MC cluster found
    if(par.mc) {
      // Find top 2 sim particles by total energy deposited
      std::map<const SimParticle*, float> sim_energy_map;
      const auto& edeps = par.mc->energyDeposits();
      for(const auto& edep : edeps) {
        if(edep.sim().isNonnull()) {
          sim_energy_map[&(*edep.sim())] += edep.energyDep();
        }
      }
      // Find primary sim
      const SimParticle* top1 = nullptr;
      float top1_energy = 0.f;
      for(const auto& [sim, energy] : sim_energy_map) {
        if(energy > top1_energy) {
          top1 = sim;
          top1_energy = energy;
        }
      }

      // Find the second sim with the highest energy that is not related to the primary sim
      const SimParticle* top2 = nullptr;
      float top2_energy = 0.f;
      for(const auto& [sim, energy] : sim_energy_map) {
        if(isRelated(sim, top1) || isRelated(top1, sim)) continue;
        if(energy > top2_energy) {
          top2 = sim;
          top2_energy = energy;
        }
      }
      par.primary_sim = top1;
      par.secondary_sim = top2;
    }
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::initLinePar(LinePar_t& par, const KalSeed* line) {
    par.init(line);
    if(!line) return;

    if(line_seed_assn_) {
      for(const auto& ent : *line_seed_assn_) {
        const art::Ptr<KalSeed>& linePtr = ent.first;
        const art::Ptr<CosmicTrackSeed>& seedPtr = ent.second;
        if(linePtr.isNonnull() && seedPtr.isNonnull()) {
          if(&(*line) == &(*linePtr)) {
            par.cosmic_seed = &(*seedPtr);
            break;
          }
        }
      }
    }
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::initCosmicSeedPar(CosmicSeedPar_t& par, const CosmicTrackSeed* seed) {
    par.init(seed);
    if(!seed) return;
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::initTimeClusterPar(TimeClusterPar_t& par, const TimeCluster* tc) {
    par.init(tc);
    if(!tc) return;

    if(sim_par_.sim && mc_digi_col_) {
      // Count the number of primary digis in the time cluster
      int n_primary_digis = 0;
      int n_other_digis = 0;
      const SimParticle& sim = *sim_par_.sim;
      std::vector<StrawHitIndex> shiv;
      const auto hit_indices = tc->hits();
      for(size_t i_hit = 0; i_hit < hit_indices.size(); ++i_hit) {
        const size_t hit_index = hit_indices.at(i_hit);
        if(hit_index >= combo_hit_col_->size()) {
          std::cerr << "[Run1BAna::" << __func__ << "] Warning: hit index " << hit_index << " out of range for combo hit collection of size "
                    << combo_hit_col_->size()
                    << " (i_hit = " << i_hit << ", N(hit indices) = " << hit_indices.size() << ")"
                    << std::endl;
          continue;
        }
        shiv.clear();
        combo_hit_col_->fillStrawHitIndices(hit_index, shiv);
        for(const size_t digi_index : shiv) {
          if(digi_index >= mc_digi_col_->size()) {
            std::cerr << "[Run1BAna::" << __func__ << "] Warning: digi index " << digi_index << " out of range for MC digi collection of size " << mc_digi_col_->size() << std::endl;
            continue;
          }
          const auto& digi = mc_digi_col_->at(digi_index);
          if(!digi.containsSimulation()) continue;
          const auto sim_ptr = digi.strawGasStep(digi.earlyEnd())->simParticle();
          if(sim_ptr.isNull()) continue;
          if(&(*sim_ptr) == &sim) {
            ++n_primary_digis;
          } else {
            ++n_other_digis;
          }
        }
      }
      par.n_primary_hits = n_primary_digis;
      par.n_other_hits = n_other_digis;
      par.n_total_primary_hits = sim_par_.nhits;
    }
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::matchLineToCluster(ClusterPar_t& par, const KalSeedCollection* lines) {
    if(!par.cluster) return;
    par.line = nullptr;
    if(!lines) return;
    const auto cluster = par.cluster;

    constexpr double max_dt = 100.; // ns
    constexpr double max_dr = 100.; // mm
    float dt_curr(1.e10), dr_curr(1.e10);
    for(const auto& line : *lines) {
      if(line.hasCaloCluster() && &(*line.caloCluster()) == &(*cluster)) ++par.nfit_matched_lines; // the fit connected the cluster and line
      if(!isGoodLine(&line)) continue;

      // check for its agreement with the cluster
      // evaluate the line distance at the calorimeter
      const auto line_pos_t = lineAtCluster(cluster, &line);
      const auto line_pos = line_pos_t.vect();
      const auto line_t = line_pos_t.t();
      const auto cl_pos = cluster->cog3Vector();
      const double dx = line_pos.x() - cl_pos.x();
      const double dy = line_pos.y() - cl_pos.y();
      const double dr = std::sqrt(dx*dx + dy*dy);
      const float dt = std::abs(cluster->time() - line_t);

      if(dr > max_dr || dt > max_dt) {
        if(debug_level_ > 0 &&
           (dr < max_dr || dt < max_dt) &&
           cluster->energyDep() > 70.) {
          std::cout << "[Run1BAna::" << __func__ << "] " << event_->id()
                    << ": Line matched in one of space / time:" << std::endl
                    << "  dt = " << dt << " dr = " << dr << std::endl;
        }
        continue;
      }
      ++par.nmatched_lines; // a line was matched to the cluster
      if(dt < dt_curr) {
        par.line = &line;
        dt_curr = dt;
        dr_curr = dr;
      }
    }
    if(par.line && debug_level_ > 2) std::cout << "[Run1BAna::" << __func__ << "] "
                                               << " Matched line with dr = " << dr_curr << " and dt = " << dt_curr
                                               << std::endl;
    if(par.line && par.nmatched_lines <= 0) std::cout << "[Run1BAna::" << __func__ << "] "
                                                      << " Matched line but N(matched lines) <= 0 = " << par.nmatched_lines
                                                      << std::endl;
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::matchSeedToCluster(ClusterPar_t& par, const CosmicTrackSeedCollection* seeds) {
    if(!par.cluster) return;
    par.cosmic_seed = nullptr;
    if(!seeds) return;
    const auto cluster = par.cluster;

    CLHEP::Hep3Vector cl_pos(cluster->cog3Vector());

    constexpr double max_dt = 100.; // ns
    constexpr double max_dr = 500.; // mm
    float dt_curr(1.e10), dr_curr(1.e10);
    for(const auto& seed : *seeds) {
      if(seed.hasCaloCluster() && &(*seed.caloCluster()) == &(*cluster)) ++par.nfit_matched_cosmic_seeds; // the fit connected the cluster and seed
      if(!isGoodCosmicSeed(&seed)) continue;

      const auto seed_pos_t = lineSeedAtCluster(cluster, &seed);
      const auto seed_pos = seed_pos_t.vect();
      const auto seed_t = seed_pos_t.t();
      const double dx = seed_pos.x() - cl_pos.x();
      const double dy = seed_pos.y() - cl_pos.y();
      const double dr = std::sqrt(dx*dx + dy*dy);
      if(dr > max_dr) {
        continue;
      }

      // check for its agreement with the cluster
      const double dt = std::abs(cluster->time() - seed_t);
      if(dt > max_dt) {
        continue;
      }
      ++par.nmatched_cosmic_seeds; // a cosmic seed was matched to the cluster
      if(dt < dt_curr) {
        par.cosmic_seed = &seed;
        dt_curr = dt;
        dr_curr = dr;
      }
    }
    if(par.cosmic_seed && debug_level_ > 2) std::cout << "[Run1BAna::" << __func__ << "] "
                                                      << " Matched cosmic seed with dt = " << dt_curr
                                                      << " and dr = " << dr_curr
                                                      << std::endl;
  }
  //--------------------------------------------------------------------------------------
  void Run1BAna::matchTimeClusterToCluster(ClusterPar_t& par, const TimeClusterCollection* time_clusters) {
    if(!par.cluster) return;
    par.time_cluster = nullptr;
    if(!time_clusters) return;
    const auto cluster = par.cluster;

    constexpr double max_dt = 100.; // ns
    constexpr double max_dr = 1000.; // mm
    float dt_curr(1.e10), dr_curr(1.e10);
    for(const auto& time_cluster : *time_clusters) {
      if(time_cluster.hasCaloCluster() && &(*time_cluster.caloCluster()) == &(*cluster)) ++par.nfit_matched_time_clusters; // the fit connected the cluster and time_cluster
      if(debug_level_ > 2) std::cout << "[Run1BAna::" << __func__ << "] "
                                     << " Checking time cluster with t0 = " << time_cluster.t0().t0() << " ns and position = " << time_cluster.position()
                                     << " N(fit matched clusters) = " << par.nfit_matched_time_clusters
                                     << std::endl;
      if(!isGoodTimeCluster(&time_cluster)) continue;

      // check for its agreement with the cluster
      const auto tc_pos = time_cluster.position();
      const auto tc_t = time_cluster.t0().t0();
      const auto cl_pos = cluster->cog3Vector();
      const double dx = tc_pos.x() - cl_pos.x();
      const double dy = tc_pos.y() - cl_pos.y();
      const double dr = std::sqrt(dx*dx + dy*dy);
      const float dt = std::abs(cluster->time() - tc_t);

      if(dr > max_dr || dt > max_dt) {
        continue;
      }
      ++par.nmatched_time_clusters; // a time cluster was matched to the cluster
      if(dt < dt_curr) {
        if(debug_level_ > 2) std::cout << "[Run1BAna::" << __func__ << "] "
                                       << " Found better time cluster match with dr = " << dr << " and dt = " << dt
                                       << std::endl;
        par.time_cluster = &time_cluster;
        dt_curr = dt;
        dr_curr = dr;
      }
    }
    if(par.time_cluster && debug_level_ > 2) std::cout << "[Run1BAna::" << __func__ << "] "
                                                       << " Matched time cluster with dr = " << dr_curr << " and dt = " << dt_curr
                                                       << std::endl;
  }

  //--------------------------------------------------------------------------------------
  float Run1BAna::getTotalEnergyDepositedBySim(const CaloClusterMC* mc, const SimParticle* sim) {
    if(!mc || !sim) return 0.f;
    float total_edep = 0.f;
    const auto& edeps = mc->energyDeposits();
    for(const auto& edep : edeps) {
      //if(edep.sim().isNonnull() && &(*edep.sim()) == sim) {
      if(edep.sim().isNonnull() && isRelated(&(*edep.sim()), sim)) { // include related sims in energy sum
        total_edep += edep.energyDep();
      }
    }
    return total_edep;
  }

  // Return the energy-weighted average time of deposits for a given sim in a CaloClusterMC
  float Run1BAna::getAverageTimeDepositedBySim(const CaloClusterMC* mc, const SimParticle* sim) {
    if(!mc || !sim) return 0.f;
    double weighted_time = 0.;
    double total_edep = 0.;
    const auto& edeps = mc->energyDeposits();
    for(const auto& edep : edeps) {
      if(edep.sim().isNonnull() && isRelated(&(*edep.sim()), sim)) { // include related sims in time average
        const double e = static_cast<double>(edep.energyDep());
        weighted_time += e * static_cast<double>(edep.time());
        total_edep += e;
      }
    }
    if(total_edep > 0.) return static_cast<float>(weighted_time / total_edep);
    return 0.f;
  }

  // Return the energy-weighted average time of all deposits in a CaloClusterMC
  float Run1BAna::getAverageTimeDeposited(const CaloClusterMC* mc) {
    if(!mc) return 0.f;
    double weighted_time = 0.;
    double total_edep = 0.;
    const auto& edeps = mc->energyDeposits();
    for(const auto& edep : edeps) {
      if(edep.sim().isNonnull()) {
        const double e = static_cast<double>(edep.energyDep());
        weighted_time += e * static_cast<double>(edep.time());
        total_edep += e;
      }
    }
    if(total_edep > 0.) return static_cast<float>(weighted_time / total_edep);
    return 0.f;
  }

  // Compute the total energy and energy-weighted average time for deposits
  // in a CaloShowerSimCollection associated to a given SimParticle (includes related sims)
  void Run1BAna::getShowerSimEnergyAndAvgTime(const CaloShowerSimCollection* col, const SimParticle* sim, float& energy, float& avgTime) const {
    energy = 0.f;
    avgTime = 0.f;
    if(!col || !sim) return;
    double tnum = 0.;
    double tsum = 0.;
    for(const auto& s : *col) {
      if(!s.sim().isNonnull()) continue;
      const SimParticle* ssim = &(*s.sim());
      if(isRelated(ssim, sim) || isRelated(sim, ssim)) {
        const double e = static_cast<double>(s.energyDep());
        tnum += e * static_cast<double>(s.time());
        tsum += e;
      }
    }
    if(tsum > 0.) {
      energy = static_cast<float>(tsum);
      avgTime = static_cast<float>(tnum / tsum);
    }
  }

  CLHEP::Hep3Vector Run1BAna::getCrystalPosition(const int crystalID) const {
    mu2e::GeomHandle<mu2e::Calorimeter> cal;
    CLHEP::Hep3Vector pos = cal->crystal(crystalID).position();
    cal->geomUtil().crystalToMu2e(crystalID, pos);
    return pos;
  }

  // Compute energy-weighted average position of crystal hits for a given SimParticle
  CLHEP::Hep3Vector Run1BAna::getSimParticleHitPosition(const CaloShowerSimCollection* col, const SimParticle* sim) const {
    CLHEP::Hep3Vector weightedPos(0., 0., 0.);
    double totalEnergy = 0.;
    if(!col || !sim) return weightedPos;

    for(const auto& shower : *col) {
      if(shower.sim().isNull()) continue;
      const auto shower_sim = &(*shower.sim());
      if(isRelated(shower_sim, &(*sim)) || isRelated(&(*sim), shower_sim)) {
        const float energy = shower.energyDep();
        const int crystalID = shower.crystalID();
        const CLHEP::Hep3Vector crystalPos = getCrystalPosition(crystalID);
        if(debug_level_ > 3) std::cout << __func__ << ": Hit position for crystal " << crystalID << ": " << crystalPos
                                       << std::endl;
        weightedPos += energy * crystalPos;
        totalEnergy += energy;
      }
    }

    if(totalEnergy > 0.) {
      weightedPos *= (1.0 / totalEnergy);
      if(debug_level_ > 3) std::cout << __func__ << ": --> Weighted position: " << weightedPos << std::endl;
    } else {
      weightedPos.set(-1000., -1000., -1000.); // make clear no hits were associated to this sim particle
      if(debug_level_ > 3) std::cout << __func__ << ": --> Default position: " << weightedPos << std::endl;
    }
    return weightedPos;
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::analyze(const art::Event& event){
    watch_->Increment("Total");
    watch_->SetTime("Event");
    ++nevt_;
    hist_norm_->Fill(0.);
    event_ = &event;
    cluster_par_.calorimeter = calorimeter_;

    //--------------------------------------------------------------------------------------
    // Retrieve the collections
    //--------------------------------------------------------------------------------------

    watch_->SetTime("DataRetrieval");
    auto clusterH = event.getValidHandle<CaloClusterCollection>(clusters_tag_); // require clusters and a trigger
    auto triggerH = event.getValidHandle<art::TriggerResults>  (trig_tag_);
    art::Handle<SimParticleCollection> simH                 ; event.getByLabel(sim_tag_            , simH);
    art::Handle<PrimaryParticle>       primaryH             ; event.getByLabel(primary_tag_        , primaryH);
    art::Handle<StrawDigiMCCollection> mc_digiH             ; event.getByLabel(mc_digi_tag_        , mc_digiH);
    art::Handle<ComboHitCollection>    combo_hitsH          ; event.getByLabel(combo_hits_tag_     , combo_hitsH);
    art::Handle<CaloHitCollection>     calo_hitsH           ; event.getByLabel(calo_hits_tag_      , calo_hitsH);
    art::Handle<CaloClusterMCCollection> calo_cluster_mcH   ; event.getByLabel(calo_cluster_mc_tag_, calo_cluster_mcH);
    art::Handle<CaloClusterMCTruthAssn> calo_cluster_mcassnH; event.getByLabel(calo_cluster_mc_tag_, calo_cluster_mcassnH);
    art::Handle<CaloShowerSimCollection> calo_shower_simH   ; event.getByLabel(calo_shower_sim_tag_, calo_shower_simH);
    art::Handle<KalSeedCollection>     lineH                ; event.getByLabel(line_tag_           , lineH);
    art::Handle<CosmicTrackSeedCollection> cosmic_seedH     ; event.getByLabel(cosmic_seed_tag_    , cosmic_seedH);
    art::Handle<KalLineAssns>              line_seed_assnH  ; event.getByLabel(line_tag_           , line_seed_assnH);
    art::Handle<TimeClusterCollection> time_clusterH        ; event.getByLabel(time_cluster_tag_   , time_clusterH);
    art::Handle<ProtonBunchIntensity>  pbiH                 ; event.getByLabel(pbi_tag_            , pbiH);

    TriggerResultsNavigator trigNav(triggerH.product());
    sim_col_              = (simH                .isValid()) ? simH.product()                 : nullptr;
    primary_              = (primaryH            .isValid()) ? primaryH.product()             : nullptr;
    mc_digi_col_          = (mc_digiH            .isValid()) ? mc_digiH.product()             : nullptr;
    combo_hit_col_        = (combo_hitsH         .isValid()) ? combo_hitsH.product()          : nullptr;
    calo_hit_col_         = (calo_hitsH          .isValid()) ? calo_hitsH.product()           : nullptr;
    calo_cluster_mc_col_  = (calo_cluster_mcH    .isValid()) ? calo_cluster_mcH.product()     : nullptr;
    calo_cluster_mc_assn_ = (calo_cluster_mcassnH.isValid()) ? calo_cluster_mcassnH.product() : nullptr;
    calo_shower_sim_col_  = (calo_shower_simH    .isValid()) ? calo_shower_simH.product()     : nullptr;
    line_col_             = (lineH               .isValid()) ? lineH.product()                : nullptr;
    cosmic_seed_col_      = (cosmic_seedH        .isValid()) ? cosmic_seedH.product()         : nullptr;
    line_seed_assn_       = (line_seed_assnH     .isValid()) ? line_seed_assnH.product()      : nullptr;
    time_cluster_col_     = (time_clusterH       .isValid()) ? time_clusterH.product()        : nullptr;
    cluster_col_          = clusterH.product();
    trig_nav_             = &trigNav;
    watch_->StopTime("DataRetrieval");

    if(debug_level_ > 1) std::cout << "[Run1BAna::" << __func__ << "::" << moduleDescription().moduleLabel() << "]"
                                   << " Input from:"
                                   << "\n  " << clusters_tag_.encode().c_str()   << ": N(clusters) = "   << cluster_col_->size()
                                   << "\n  " << calo_hits_tag_.encode().c_str()  << ": N(calo hits) = "  << int((calo_hit_col_)   ? calo_hit_col_->size()  : -1)
                                   << "\n  " << combo_hits_tag_.encode().c_str() << ": N(combo hits) = " << int((combo_hit_col_)  ? combo_hit_col_->size() : -1)
                                   << "\n  " << mc_digi_tag_.encode().c_str()    << ": N(MC digis) = "   << int((mc_digi_col_)    ? mc_digi_col_->size()    : -1)
                                   << "\n  " << calo_cluster_mc_tag_.encode().c_str() << ": N(calo cluster MC) = " << int((calo_cluster_mc_col_) ? calo_cluster_mc_col_->size() : -1)
                                   << "\n  " << time_cluster_tag_.encode().c_str() << ": N(time clusters) = " << int((time_cluster_col_) ? time_cluster_col_->size() : -1)
                                   << "\n  " << line_tag_.encode().c_str()       << ": N(lines) = "      << int((line_col_)       ? line_col_->size()       : -1)
                                   << std::endl;

    if(debug_level_ > 2 && time_cluster_col_) {
      std::cout << "[Run1BAna::" << __func__ << "::" << moduleDescription().moduleLabel() << "]"
                << " Time clusters:" << std::endl;
      for(const auto& tc : *time_cluster_col_) {
        std::cout << "  t0 = " << tc.t0().t0() << " ns, position = " << tc.position() << ", N(hits) = " << tc.nhits() << std::endl;
        for(const auto& hit_index : tc.hits()) {
          std::cout << "    hit index = " << hit_index << std::endl;
        }
      }
    }
    //--------------------------------------------------------------------------------------
    // Event filtering
    //--------------------------------------------------------------------------------------

    if(max_gen_energy_ > 0.) {
      if(!primary_) throw std::runtime_error("Requiring a gen energy cut but no primary is found!");
      const SimParticle* primary_sim = (!primary_->primarySimParticles().empty()) ? &(*primary_->primarySimParticles().front()) : nullptr;
      if(!primary_sim) throw std::runtime_error("Requiring a gen energy cut but no primary sim is found!");
      const double gen_energy = primary_sim->startMomentum().e();
      if(gen_energy > max_gen_energy_) {
        watch_->StopTime("Event");
        return;
      }
    }

    //--------------------------------------------------------------------------------------
    // Initialize fields
    //--------------------------------------------------------------------------------------

    SimUtils::fillSimInfo(sim_info_, sim_col_, mc_digi_col_);
    cluster_par_.init();
    line_par_.init();
    cosmic_seed_par_.init();
    time_cluster_par_.init();
    photon_.init();
    evt_par_.init((pbiH.isValid()) ? pbiH->intensity() : 0);
    const SimParticle* primary_sim = (primary_ && !primary_->primarySimParticles().empty()) ? &(*primary_->primarySimParticles().front()) : nullptr;
    const SimParticle* secondary_sim = (primary_ && primary_->primarySimParticles().size() > 1) ? &(*primary_->primarySimParticles().at(1)) : nullptr;
    sim_par_.init((primary_sim) ? primary_sim : nullptr, (primary_sim) ? sim_info_[primary_sim->id().asUint()].nhits_ : 0, 0.);

    // If the first event, initialize the trigger path histogram bins for stability
    if(nevt_ == 1) {
      for(int iset = 0; iset < kMaxHists; ++iset) {
        if(!hist_[iset]) continue;
        TH1* h = hist_[iset]->event.trig_paths;
        if(!h) continue;
        for (size_t index = 0; index< trigNav.getTrigPaths().size(); ++index) {
          const std::string path = trigNav.getTrigPathNameByIndex(index);
          h->GetXaxis()->SetBinLabel(h->FindBin(index), path.c_str());
        }
      }
    }

    //--------------------------------------------------------------------------------------
    // First loop over collections to count objects with IDs, find global objects
    //--------------------------------------------------------------------------------------

    watch_->SetTime("CountObjects");
    // Clusters
    const CaloCluster* max_cluster = nullptr;
    const CaloCluster* best_cluster = nullptr;
    for(const auto& cluster : *(cluster_col_)) {
      // initClusterPar(cluster_par_, &cluster);
      // line_par_.init(cluster_par_.line);

      const bool cluster_time = (cluster.time() > 600. && cluster.time() < 1650.);
      const bool cluster_id = isGoodCluster(&cluster) && (cluster.energyDep() > 70. && cluster_time);

      if(cluster_id) {
        // Count good clusters
        ++evt_par_.n_good_clusters;

        // Highest energy "good" cluster
        if(!best_cluster || best_cluster->energyDep() < cluster.energyDep()) best_cluster = &cluster;
      }

      // Find the highest energy cluster
      if(!max_cluster || max_cluster->energyDep() < cluster.energyDep()) max_cluster = &cluster;

    }

    // Lines
    if(line_col_) {
      for(const auto& line : *(line_col_)) {
        // line_par_.init(&line);
        // // initialize the associated cluster if defined
        // const CaloCluster* cluster = (line.hasCaloCluster()) ? &(*line.caloCluster()) : nullptr;
        // initClusterPar(cluster_par_, cluster);
        // cluster_par_.line = &line;

        bool line_id = isGoodLine(&line);
        if(line_id) ++evt_par_.n_good_lines;
      }
    }

    // Cosmic track seeds
    if(cosmic_seed_col_) {
      for(const auto& seed : *cosmic_seed_col_) {
        // cosmic_seed_par_.init(&seed);
        if(isGoodCosmicSeed(&seed)) {
          ++evt_par_.n_good_cosmic_seeds;
        }
      }
    }

    // Time clusters
    if(time_cluster_col_) {
      for(const auto& time_cluster : *time_cluster_col_) {
        // initTimeClusterPar(time_cluster_par_,&time_cluster);
        time_cluster_par_.init(&time_cluster);
        if(isGoodTimeCluster(&time_cluster)) {
          ++evt_par_.n_good_time_clusters;
        }
      }
    }
    watch_->StopTime("CountObjects");

    //--------------------------------------------------------------------------------------
    // Fill histograms
    //--------------------------------------------------------------------------------------

    watch_->SetTime("Analysis-Clusters");
    for(const auto& cluster : *(cluster_col_)) {
      watch_->SetTime("Analysis-Clusters-init");
      initClusterPar(cluster_par_, &cluster);
      initLinePar(line_par_, cluster_par_.line);
      initCosmicSeedPar(cosmic_seed_par_, (line_par_.cosmic_seed) ? line_par_.cosmic_seed : cluster_par_.cosmic_seed);
      initTimeClusterPar(time_cluster_par_, cluster_par_.time_cluster);
      watch_->StopTime("Analysis-Clusters-init");

      // Find the "best" cluster in the event
      bool cluster_time = (cluster.time() > 600. && cluster.time() < 1650.);
      bool cluster_id = (cluster.energyDep() > 70. && cluster_time);

      // Two primary event merged into one cluster
      bool merged_cluster = (cluster_par_.primary_sim && cluster_par_.secondary_sim &&
                             (isRelated(cluster_par_.primary_sim, primary_sim) || isRelated(cluster_par_.primary_sim, secondary_sim)) &&
                             (isRelated(cluster_par_.secondary_sim, primary_sim) || isRelated(cluster_par_.secondary_sim, secondary_sim)));

      // All clusters
      fillHistograms(hist_[1]);

      // Clusters above 70 MeV
      if(cluster.energyDep() > 70.) {
        fillHistograms(hist_[2]);
        if(cluster.diskID() == 0) fillHistograms(hist_[3]);
        else                      fillHistograms(hist_[4]);
      }
      if(cluster_id) fillHistograms(hist_[5]);

      // Clusters above thresholds and in time
      if(cluster_time) {
        if(cluster.energyDep() > 30.) fillHistograms(hist_[6]);
        if(cluster.energyDep() > 50.) fillHistograms(hist_[7]);
        if(merged_cluster) {
          fillHistograms(hist_[30]);
          const float edep_1 = getTotalEnergyDepositedBySim(cluster_par_.mc, cluster_par_.primary_sim);
          const float edep_2 = getTotalEnergyDepositedBySim(cluster_par_.mc, cluster_par_.secondary_sim);
          if(edep_1 > 10. && edep_2 > 10.) {
            std::cout << "Merged cluster: " << event.id() << std::endl;
            fillHistograms(hist_[31]);
          }
        } else {
          fillHistograms(hist_[40]);
          if(cluster.energyDep() > 30.) {
            fillHistograms(hist_[41]);
            if(primary_ && calo_shower_sim_col_ && primary_->primarySimParticles().size() > 1) {
              const auto& pvec = primary_->primarySimParticles();
              const SimParticle* psim0 = &(*pvec.front());
              const SimParticle* psim1 = &(*pvec.at(1));
              float e0 = 0.f, t0 = 0.f;
              float e1 = 0.f, t1 = 0.f;
              getShowerSimEnergyAndAvgTime(calo_shower_sim_col_, psim0, e0, t0);
              getShowerSimEnergyAndAvgTime(calo_shower_sim_col_, psim1, e1, t1);
              const auto pos0 = getSimParticleHitPosition(calo_shower_sim_col_, psim0);
              const auto pos1 = getSimParticleHitPosition(calo_shower_sim_col_, psim1);
              const double dx = pos0.x() - pos1.x();
              const double dy = pos0.y() - pos1.y();
              const double dz = pos0.z() - pos1.z();
              const double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
              if(dist < 50. && std::fabs(t0-t1) < 25.)
                std::cout << "Unmerged cluster, dr = " << dist << " dt = " << t0 - t1 << ": " << event.id() << std::endl;
            }
          }
        }
      }
    }
    watch_->StopTime("Analysis-Clusters");

    // Per-line (KalSeed) histograms
    watch_->SetTime("Analysis-Lines");
    if(line_col_) {
      for(const auto& line : *(line_col_)) {
        initLinePar(line_par_, &line);
        initCosmicSeedPar(cosmic_seed_par_, line_par_.cosmic_seed);
        initTimeClusterPar(time_cluster_par_, line_par_.time_cluster);
        // initialize the associated cluster if defined
        const CaloCluster* cluster = (line.hasCaloCluster()) ? &(*line.caloCluster()) : nullptr;
        initClusterPar(cluster_par_, cluster);
        cluster_par_.line = &line;

        // fill all histograms
        fillHistograms(hist_[80]);
        if(cluster) {
          fillHistograms(hist_[81]);
          const bool id = cluster->energyDep() > 70. && cluster->time() < 1650. && cluster->time() > 600.;
          if(id) fillHistograms(hist_[82]);
        }
      }
    }
    watch_->StopTime("Analysis-Lines");

    // Per cosmic seed histograms
    watch_->SetTime("Analysis-CosmicSeeds");
    if(cosmic_seed_col_) {
      for(const auto& seed : *cosmic_seed_col_) {
        initCosmicSeedPar(cosmic_seed_par_, &seed);
        line_par_.init(nullptr); // for now, no line association to cosmic seeds
        time_cluster_par_.init(nullptr); // for now, no time cluster association to cosmic seeds
        if(seed.hasCaloCluster()) initClusterPar(cluster_par_, &(*seed.caloCluster()));
        else                      initClusterPar(cluster_par_, nullptr);
        fillHistograms(hist_[90]);
        if(isGoodCosmicSeed(&seed)) {
          fillHistograms(hist_[91]);
          if(cluster_par_.cluster) {
            fillHistograms(hist_[92]);
            if(cluster_par_.cluster->energyDep() > 70.) {
              fillHistograms(hist_[93]);
            }
          }
        }
      }
    }
    watch_->StopTime("Analysis-CosmicSeeds");

    // Per time cluster histograms
    watch_->SetTime("Analysis-TimeClusters");
    if(time_cluster_col_) {
      const TimeCluster* max_sim_hits_tc = nullptr;
      int max_sim_hits = -1;
      for(const auto& time_cluster : *time_cluster_col_) {
        initTimeClusterPar(time_cluster_par_, &time_cluster);
        line_par_.init(nullptr); // for now, no line association to time clusters
        cosmic_seed_par_.init(nullptr); // for now, no cosmic seed association to time clusters
        if(time_cluster.hasCaloCluster()) initClusterPar(cluster_par_, &(*time_cluster.caloCluster()));
        else                              initClusterPar(cluster_par_, nullptr);
        fillHistograms(hist_[95]);
        if(time_cluster_par_.n_primary_hits > max_sim_hits) {
          max_sim_hits = time_cluster_par_.n_primary_hits;
          max_sim_hits_tc = &time_cluster;
        }
      }
      if(max_sim_hits_tc && max_sim_hits > 0) { // must at least be 1 hit
        initTimeClusterPar(time_cluster_par_, max_sim_hits_tc);
        line_par_.init(nullptr); // for now, no line association to time clusters
        cosmic_seed_par_.init(nullptr); // for now, no cosmic seed association to time clusters
        if(max_sim_hits_tc->hasCaloCluster()) initClusterPar(cluster_par_, &(*max_sim_hits_tc->caloCluster()));
        else                                  initClusterPar(cluster_par_, nullptr);
        fillHistograms(hist_[96]);
      }
    }

    // All events, highest energy cluster
    watch_->SetTime("Analysis-Event");
    initClusterPar(cluster_par_, max_cluster);
    initLinePar(line_par_, cluster_par_.line);
    initCosmicSeedPar(cosmic_seed_par_, (line_par_.cosmic_seed) ? line_par_.cosmic_seed : cluster_par_.cosmic_seed);
    initTimeClusterPar(time_cluster_par_, cluster_par_.time_cluster);
    fillHistograms(hist_[0]);

    // above 70 MeV
    if(cluster_par_.cluster && cluster_par_.cluster->energyDep() > 70.) {
      fillHistograms(hist_[10]);
      const float dt = cluster_par_.cluster->time() - lineAtCluster(cluster_par_.cluster, cluster_par_.line).t();
      if(std::fabs(dt) < 50.) fillHistograms(hist_[11]);
      else                    fillHistograms(hist_[12]);
    }

    // at least 10 sim hits
    if(sim_par_.nhits >= 10) {
      fillHistograms(hist_[8]);
      if(cluster_par_.cluster && cluster_par_.cluster->energyDep() > 70.) fillHistograms(hist_[9]);
    }

    // "best" cluster in the event, passing the cluster ID
    if(best_cluster) {
      initClusterPar(cluster_par_, best_cluster);
      initLinePar(line_par_, cluster_par_.line);
      initCosmicSeedPar(cosmic_seed_par_, (line_par_.cosmic_seed) ? line_par_.cosmic_seed : cluster_par_.cosmic_seed);
      initTimeClusterPar(time_cluster_par_, cluster_par_.time_cluster);
      fillHistograms(hist_[20]);
      const float dt = cluster_par_.cluster->time() - lineAtCluster(cluster_par_.cluster, cluster_par_.line).t();
      if(std::fabs(dt) < 50.) fillHistograms(hist_[21]);
      else                    fillHistograms(hist_[22]);
      if(best_cluster->energyDep() > 80.) fillHistograms(hist_[23]);

      // Signal selection
      bool signal_id = true;
      signal_id &= cluster_par_.ncr() > 0;
      signal_id &= cluster_par_.ncr() < 6;
      signal_id &= cluster_par_.frac_1() > 0.60;
      signal_id &= cluster_par_.frac_2() > 0.80;
      signal_id &= cluster_par_.t_var < 1.0;
      signal_id &= cluster_par_.second_moment < 1.e5;
      // if(cluster_par_.time_cluster) signal_id &= cluster_par_.time_cluster->nhits() < 80; // FIXME: Depends on the beam intensity used
      if(signal_id) fillHistograms(hist_[25]);
    }

    // histograms for photon candidates
    photon_.init(best_cluster); // FIXME: Find all photon candidates
    if(photon_.cluster) {
      const double prev_weight = evt_par_.weight; // store the previous weight
      initClusterPar(cluster_par_, photon_.cluster);
      initLinePar(line_par_, cluster_par_.line);
      initCosmicSeedPar(cosmic_seed_par_, (line_par_.cosmic_seed) ? line_par_.cosmic_seed : cluster_par_.cosmic_seed);
      initTimeClusterPar(time_cluster_par_, cluster_par_.time_cluster);
      photon_.cosmic_seed = cluster_par_.cosmic_seed;
      photon_.line = cluster_par_.line;
      photon_.time_cluster = cluster_par_.time_cluster;
      const float radius = photon_.r;
      const auto cluster = photon_.cluster;
      const int disk_id = cluster->diskID();

      // Get the RMC weight, if this is a flat photon sample
      double weight = 1.;
      if(primary_sim && primary_sim->creationCode() == mu2e::ProcessCode::mu2eFlatPhoton) {
        const double e_gen = primary_sim->startMomentum().e();
        constexpr double kmax = 90.1; // use for reference
        weight = Run1BAnaUtils::closureApprox(e_gen, kmax);
      }

      // Assign the event weight
      evt_par_.weight = weight;
      fillHistograms(hist_[60]);

      // Set without weights
      evt_par_.weight = 1.;
      fillHistograms(hist_[61]);
      evt_par_.weight = weight;

      if(!photon_.time_cluster) { // no time cluster selection
        fillHistograms(hist_[62]);
        if(disk_id == 0 && radius > 500.) fillHistograms(hist_[67]);
        if(disk_id == 0 && radius > 550.) fillHistograms(hist_[68]);
      } else {
        const float dt = photon_.cluster->time() - photon_.time_cluster->t0().t0();
        if(std::fabs(dt) < 50.) fillHistograms(hist_[63]);
        else                    fillHistograms(hist_[64]);
      }
      // line/seed selections
      if(!photon_.cosmic_seed) fillHistograms(hist_[65]);
      if(!photon_.line)        fillHistograms(hist_[66]);

      evt_par_.weight = prev_weight;
    }

    watch_->StopTime("Analysis-Event");


    watch_->StopTime("Event");
  } // end analyze


  //--------------------------------------------------------------------------------------
  void Run1BAna::endRun(const art::Run& run) {
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::endJob() {
    watch_->StopTime("Total");
    watch_->Print(std::cout);
  }

}
using mu2e::Run1BAna;
DEFINE_ART_MODULE(Run1BAna)
