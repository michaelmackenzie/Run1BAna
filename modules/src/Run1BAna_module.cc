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
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloClusterMC.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/Mu2eUtilities/inc/StopWatch.hh"
#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"
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
#include "Run1BAna/modules/inc/SimUtils.hh"

using namespace CLHEP;

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
      fhicl::Atom<art::InputTag>      simCol          { Name("simCollection")          , Comment("simCollection")                       };
      fhicl::Atom<art::InputTag>      primary         { Name("primaryParticle")        , Comment("primaryParticle")                     };
      fhicl::Atom<art::InputTag>      digiCol         { Name("strawDigiMCCollection")  , Comment("strawDigiMCCollection")               };
      fhicl::Atom<art::InputTag>      comboCol        { Name("ComboHitCollection")     , Comment("ComboHitCollection")                  };
      fhicl::Atom<art::InputTag>      caloHitCol      { Name("CaloHitCollection")      , Comment("CaloHitCollection")                   };
      fhicl::Atom<art::InputTag>      caloClusterCol  { Name("caloClusterCollection")  , Comment("caloClusterCollection")               };
      fhicl::Atom<art::InputTag>      caloClusterMCCol{ Name("caloClusterMCCollection"), Comment("caloClusterMCCollection")             };
      fhicl::Atom<art::InputTag>      lineCol         { Name("LineCollection")         , Comment("Kinematic line collection")           };
      fhicl::Atom<std::string>        trigProcess     { Name("triggerProcess")         , Comment("Process name for the trigger results")};
      fhicl::Atom<art::InputTag>      pbi             { Name("PBI")                    , Comment("ProtonBunchIntensity tag")            };
      fhicl::Atom<double>             maxGenEnergy    { Name("maxGenEnergy")           , Comment("Cut on the maximum primary energy")  , -1.};
      fhicl::Atom<int>                debugLevel      { Name("debugLevel")             , Comment("debugLevel")                         , 0 };
    };


    //--------------------------------------------------------------------------------------
    // Histograms
    //--------------------------------------------------------------------------------------

    // Per event
    struct EventHist_t {
      TH1* npot;
      TH1* n_mc_digis;
      TH1* ncombo_hits;
      TH1* ncalo_hits;
      TH1* nlines;
      TH1* nclusters;
      TH1* ngood_clusters;

      TH1* trig_bits;
      TH1* trig_paths;
    };

    // Per cluster
    struct ClusterHist_t {
      TH1* energy;
      TH1* mc_edep;
      TH1* time;
      TH1* radius;
      TH1* ncr;
      TH1* disk;

      TH1* line_dt;
      TH1* line_dr;

      TH1* esum;
      TH1* dt;
      TH1* dr;
      TH2* dr_vs_dt;
    };

    // Per line (KalSeed)
    struct LineHist_t {
      TH1* chi2;
      TH1* nhits;
    };

    // Per sim
    struct SimHist_t {
      TH1* pdg;
      TH1* type;
      TH1* parent_type;
      TH1* origin_type;
      TH1* nhits;
      TH1* energy_start;
      TH1* mom_cz_start;
      TH1* edep;
    };

    // Total
    struct Hist_t {
      EventHist_t   event;
      ClusterHist_t cluster;
      LineHist_t    line;
      SimHist_t     sim;
    };

    //--------------------------------------------------------------------------------------
    // Internal data structures
    //--------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------
    struct EventPar_t {
      long long npot;
      int n_good_clusters;

      EventPar_t() {
        init();
      }

      void init(long long np = 0) {
        npot = np;
        n_good_clusters = 0;
      }
    };

    //--------------------------------------------------------------------------------------
    struct SimPar_t {
      const SimParticle* sim;
      unsigned nhits;
      double edep;

      SimPar_t() {
        init();
      }

      void init(const SimParticle* s = nullptr, unsigned n = 0, double e = 0.) {
        sim = s;
        nhits = n;
        edep = e;
      }
    };

    //--------------------------------------------------------------------------------------
    struct ClusterPar_t {
      const CaloCluster* cluster;
      const KalSeed*    line;
      float r;

      ClusterPar_t() {
        init();
      }

      void init(const CaloCluster* cl = nullptr,  const KalSeed* ln = nullptr) {
        cluster = cl;
        line = ln;
        r = 0.f;
        mc = nullptr;
        if(!cl) return;

        const float x = cluster->cog3Vector().x();
        const float y = cluster->cog3Vector().y();
        r = std::sqrt(x*x + y*y);
      }

      const CaloClusterMC* mc;

      double line_dt() const {
        if(!cluster || !line) return -9999.;
        // const auto seg = line->nearestSegment(cluster_pos)
        // if(seg == line->segments().end()) return -9999.;
        // const auto cluster_pos = cluster->cog3Vector();
        // const double t0 = seg.t0Val(TrkFitFlag(TrkFitFlag::KKLine));
        double t0;
        auto seg = line->t0Segment(t0);
        if(seg == line->segments().end()) return -9999.;
        return t0 - cluster->time();
      }

      double line_dr() const {
        if(!cluster || !line) return -9999.;
        double dr(1.e10), dz(1.e10);
        const auto cluster_pos = cluster->cog3Vector();
        for(const auto& seg : line->segments()) {
          const auto pos = seg.position3();
          const double curr_dr = std::sqrt( std::pow(pos.x() - cluster_pos.x(), 2) +
                                            std::pow(pos.y() - cluster_pos.y(), 2) );
          const double curr_dz = std::abs(pos.z() - cluster_pos.z());
          if(curr_dz < dz) {
            dz = curr_dz;
            dr = curr_dr;
          }
        }
        return dr;
      }
    };

    //--------------------------------------------------------------------------------------
    struct LinePar_t {
      const KalSeed* line;

      LinePar_t() { init(); }
      void init(const KalSeed* l = nullptr) {
        line = l;
        if(!l) return;
      }
    };

    //--------------------------------------------------------------------------------------
    struct ChargedTrack_t {
      const CaloCluster* cluster;
      const KalSeed* line;

      ChargedTrack_t() { init(); }
      void init(const CaloCluster* c = nullptr, const KalSeed* l = nullptr) {
        cluster = c;
        line = l;
        if(!c || !l) return;
      }
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

    void bookClusterHistograms(ClusterHist_t* Hist, const int index, const char* name);
    void bookEventHistograms  (EventHist_t*   Hist, const int index, const char* name);
    void bookLineHistograms   (LineHist_t*    Hist, const int index, const char* name);
    void bookSimHistograms    (SimHist_t*     Hist, const int index, const char* name);
    void bookHistograms       (                     const int index, const char* name);

    void fillClusterHistograms(ClusterHist_t* Hist);
    void fillLineHistograms   (LineHist_t*    Hist);
    void fillEventHistograms  (EventHist_t*   Hist);
    void fillSimHistograms    (SimHist_t*     Hist);
    void fillHistograms       (Hist_t*        Hist);

    void initClusterPar(ClusterPar_t& par, const CaloCluster* cluster);
    void matchLineToCluster(ClusterPar_t& par, const KalSeedCollection* lines);

    //--------------------------------------------------------------------------------------
    // Data
    //--------------------------------------------------------------------------------------

    art::InputTag  sim_tag_;
    art::InputTag  primary_tag_;
    art::InputTag  mc_digi_tag_;
    art::InputTag  combo_hits_tag_;
    art::InputTag  calo_hits_tag_;
    art::InputTag  calo_cluster_mc_tag_;
    art::InputTag  clusters_tag_;
    art::InputTag  line_tag_;
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
    const ComboHitCollection*              combo_hits_col_ = nullptr;
    const CaloHitCollection*               calo_hits_col_ = nullptr;
    const CaloClusterCollection*           cluster_col_ = nullptr;
    const CaloClusterMCCollection*         calo_cluster_mc_col_ = nullptr;
    const CaloClusterMCTruthAssn*          calo_cluster_mc_assn_ = nullptr;
    const KalSeedCollection*               line_col_ = nullptr;
    const TriggerResultsNavigator*         trig_nav_ = nullptr;
    EventPar_t                             evt_par_;
    ClusterPar_t                           cluster_par_;
    LinePar_t                              line_par_;
    SimPar_t                               sim_par_;
    std::map<unsigned, SimUtils::Sim_t>    sim_info_;
    mu2e::StopWatch*                       watch_; // for timing

    unsigned  long                         nevt_;
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
    , clusters_tag_       (config().caloClusterCol())
    , line_tag_           (config().lineCol())
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

      bookHistograms(10, "max_70MeV");
      bookHistograms(11, "max_70MeV_line");
      bookHistograms(12, "max_70MeV_noline");

      bookHistograms(20, "cluster_ID");
      bookHistograms(21, "cluster_ID_line");
      bookHistograms(22, "cluster_ID_noline");

      // For timing info
      watch_ = new mu2e::StopWatch();
      watch_->Calibrate();
    }

  //--------------------------------------------------------------------------------------
  void Run1BAna::beginJob() {
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::beginRun(const art::Run & run){
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookEventHistograms(EventHist_t* Hist, const int index, const char* title) {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(std::format("evt_{}", index), title);
    Hist->npot           = dir.make<TH1F>("npot"          , "N(POT);"                     ,  100,    0.,   1.e7);
    Hist->n_mc_digis     = dir.make<TH1F>("n_mc_digis"    , "N(MC digis);"                , 1000,    0., 10000.);
    Hist->ncombo_hits    = dir.make<TH1F>("n_combo_hits"  , "N(combo hits);"              ,  500,    0., 10000.);
    Hist->ncalo_hits     = dir.make<TH1F>("n_calo_hits"   , "N(calo hits);"               ,  500,    0.,  2000.);
    Hist->nlines         = dir.make<TH1F>("n_lines"       , "N(Kinematic lines);"         ,  100,    0.,   100.);
    Hist->nclusters      = dir.make<TH1D>("nclusters"     , "N(calo clusters);"           ,  100,    0.,   100.);
    Hist->ngood_clusters = dir.make<TH1D>("ngood_clusters", "N(calo clusters | ID);"      ,  100,    0.,   100.);

    Hist->trig_bits     = dir.make<TH1D>("trig_bits"    , "Trigger bits;"               , 1000,    0.,   1000);
    Hist->trig_paths    = dir.make<TH1D>("trig_paths"   , "Trigger paths;"              ,  100,    0.,    100);

  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookClusterHistograms(ClusterHist_t* Hist, const int index, const char* title) {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(std::format("cls_{}", index), title);
    Hist->energy        = dir.make<TH1F>("energy"  , "Cluster energy;Energy (MeV);"     ,  300,    0.,    300.);
    Hist->mc_edep       = dir.make<TH1F>("mc_edep" , "Cluster MC energy;MC Energy (MeV);",  300,    0.,    300.);
    Hist->time          = dir.make<TH1F>("time"    , "Cluster time;Time (ns);"          ,  100,    0.,   2000.);
    Hist->radius        = dir.make<TH1F>("radius"  , "Cluster radial position;R (mm);"  ,  100,    0.,    700.);
    Hist->ncr           = dir.make<TH1D>("ncr"     , "Number of crystals;N(Crystals);"  ,   20,    0.,     20.);
    Hist->disk          = dir.make<TH1D>("disk"    , "Disk ID;Disk ID;"                 ,    3,   -1.,      2.);
    Hist->line_dt       = dir.make<TH1F>("line_dt" , "Line-cluster time diff;dt (ns);"  ,  100, -200.,    200.);
    Hist->line_dr       = dir.make<TH1F>("line_dr" , "Line-cluster distance;dr (mm);"   ,  100,    0.,    500.);

    Hist->esum          = dir.make<TH1F>("esum"      , "Energy sum between two clusters;E_{sum} (MeV);"    ,  300,    0.,    300.);
    Hist->dt            = dir.make<TH1F>("dt"        , "Time difference between two clusters;#Deltat (ns);",  300,    0.,    300.);
    Hist->dr            = dir.make<TH1F>("dr"        , "Position difference between two clusters;#Deltar (mm);",  300, -150.,    150.);
    Hist->dr_vs_dt      = dir.make<TH2F>("dr_vs_dt"  , "Position and time difference between two clusters;#Deltat (ns);#Deltar (mm)",  100, -100., 100., 100, 0., 200.);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookLineHistograms(LineHist_t* Hist, const int index, const char* title) {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(std::format("line_{}", index), title);
    Hist->chi2 = dir.make<TH1F>("chi2", "Line fit chi2;chi2;"               , 100, 0., 1000.);
    Hist->nhits= dir.make<TH1F>("nhits","Line nhits;N(straw hits);"         ,  50, 0.,   500.);
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
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookHistograms(const int index, const char* title) {
    if(index >= kMaxHists) throw std::runtime_error("Too many histograms!");
    hist_[index] = new Hist_t;
    auto Hist = hist_[index];
    bookEventHistograms  (&Hist->event  , index, title);
    bookClusterHistograms(&Hist->cluster, index, title);
    bookLineHistograms   (&Hist->line   , index, title);
    bookSimHistograms    (&Hist->sim    , index, title);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillClusterHistograms(ClusterHist_t* Hist) {
    if(!Hist) return;
    auto Cluster = cluster_par_.cluster;
    if(!Cluster) return;

    Hist->energy->Fill(Cluster->energyDep());
    if(cluster_par_.mc) Hist->mc_edep->Fill(cluster_par_.mc->totalEnergyDep());
    Hist->time->Fill(Cluster->time());
    Hist->radius->Fill(cluster_par_.r);
    Hist->ncr->Fill(Cluster->caloHitsPtrVector().size());
    Hist->disk->Fill(Cluster->diskID());
    Hist->line_dt->Fill(cluster_par_.line_dt());
    Hist->line_dr->Fill(cluster_par_.line_dr());

    for(const auto& cls : *cluster_col_) {
      if(&cls == &(*Cluster)) continue;
      const float dx = (Cluster->cog3Vector() - cls.cog3Vector()).x();
      const float dy = (Cluster->cog3Vector() - cls.cog3Vector()).y();
      const float dr = std::sqrt(dx*dx + dy*dy);
      const float dt = Cluster->time() - cls.time();
      Hist->esum      ->Fill(Cluster->energyDep() + cls.energyDep());
      Hist->dt        ->Fill(dt);
      Hist->dr        ->Fill(dr);
      Hist->dr_vs_dt  ->Fill(dt, dr);
    }
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillLineHistograms(LineHist_t* Hist) {
    if(!Hist) return;
    // Hist->chi2->Fill(_line_par.chi2);
    // Hist->nhits->Fill(_line_par.nhits);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillEventHistograms(EventHist_t* Hist) {
    if(!Hist) return;
    Hist->npot           ->Fill(evt_par_.npot);
    Hist->n_mc_digis     ->Fill((mc_digi_col_) ? mc_digi_col_->size() : 0);
    Hist->ncombo_hits    ->Fill((combo_hits_col_) ? combo_hits_col_->size() : 0);
    Hist->ncalo_hits     ->Fill((calo_hits_col_) ? calo_hits_col_->size() : 0);
    Hist->nlines         ->Fill((line_col_) ? line_col_->size() : 0);
    Hist->nclusters      ->Fill((cluster_col_) ? cluster_col_->size() : 0);
    Hist->ngood_clusters ->Fill(evt_par_.n_good_clusters);

    // Trigger information
    for (size_t index = 0; index < trig_nav_->getTrigPaths().size(); ++index) {
      const std::string path = trig_nav_->getTrigPathName(index);
      if(trig_nav_->accepted(path)) {
        Hist->trig_bits ->Fill(trig_nav_->findTrigPathID(path));
        Hist->trig_paths->Fill(path.c_str(), 1.);
      }
    }
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillSimHistograms(SimHist_t* Hist) {
    if(!Hist) return;
    const SimParticle* sim = sim_par_.sim;
    if(!sim) return;
    Hist->pdg         ->Fill(sim->pdgId());
    Hist->type        ->Fill(SimUtils::getSimType(sim));
    Hist->parent_type ->Fill(SimUtils::getSimType((sim->hasParent()) ? &(*sim->parent()) : nullptr));
    Hist->origin_type ->Fill(SimUtils::getSimOriginType(sim));
    Hist->energy_start->Fill(sim->startMomentum().e());
    Hist->nhits       ->Fill(sim_par_.nhits);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillHistograms(Hist_t* Hist) {
    if(!Hist) return;
    watch_->SetTime("FillHistogram");
    fillEventHistograms(&Hist->event);
    fillClusterHistograms(&Hist->cluster);
    fillLineHistograms(&Hist->line);
    fillSimHistograms(&Hist->sim);
    watch_->StopTime("FillHistogram");
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::initClusterPar(ClusterPar_t& par, const CaloCluster* cluster) {
    par.init(cluster);
    matchLineToCluster(par, line_col_);

    // Match reconstructed cluster to MC cluster using the CaloClusterMCTruthAssn
    par.mc = nullptr;
    if(calo_cluster_mc_assn_ && par.cluster) {
      for(const auto& ent : *calo_cluster_mc_assn_) {
        const art::Ptr<CaloCluster>& recoPtr = ent.first;
        const art::Ptr<CaloClusterMC>& mcPtr = ent.second;
        if(recoPtr.isNonnull() && mcPtr.isNonnull()) {
          if(par.cluster == &(*recoPtr)) {
            par.mc = &(*mcPtr);
            break;
          }
        }
      }
    }
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::matchLineToCluster(ClusterPar_t& par, const KalSeedCollection* lines) {
    if(!par.cluster) return;
    par.line = nullptr;
    if(!lines) return;
    const auto cluster = par.cluster;

    constexpr double max_dt = 200.; // ns
    constexpr double max_dr = 300.; // mm
    for(const auto& line : *lines) {
      ClusterPar_t line_par;
      line_par.init(cluster, &line);
      const float dt = std::abs(line_par.line_dt());
      if(dt > max_dt) continue;
      if(0. > max_dr) continue; // FIXME: compute dr between line and cluster
      const float dt_curr = (par.line) ? std::abs(par.line_dt()): max_dt + 1.;
      if(dt < dt_curr) par.line = &line;
    }
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::analyze(const art::Event& event){
    watch_->Increment("Total");
    watch_->SetTime("Event");
    ++nevt_;
    hist_norm_->Fill(0.);

    //--------------------------------------------------------------------------------------
    // Retrieve the collections
    //--------------------------------------------------------------------------------------

    auto clusterH = event.getValidHandle<CaloClusterCollection>(clusters_tag_); // require clusters and a trigger
    auto triggerH = event.getValidHandle<art::TriggerResults>  (trig_tag_);
    art::Handle<SimParticleCollection> simH       ; event.getByLabel(sim_tag_       , simH);
    art::Handle<PrimaryParticle>       primaryH   ; event.getByLabel(primary_tag_   , primaryH);
    art::Handle<StrawDigiMCCollection> mc_digiH   ; event.getByLabel(mc_digi_tag_   , mc_digiH);
    art::Handle<ComboHitCollection>    combo_hitsH; event.getByLabel(combo_hits_tag_, combo_hitsH);
    art::Handle<CaloHitCollection>     calo_hitsH ; event.getByLabel(calo_hits_tag_ , calo_hitsH);
    art::Handle<CaloClusterMCCollection> calo_cluster_mcH; event.getByLabel(calo_cluster_mc_tag_, calo_cluster_mcH);
    art::Handle<CaloClusterMCTruthAssn> calo_cluster_mcassnH; event.getByLabel(calo_cluster_mc_tag_, calo_cluster_mcassnH);
    art::Handle<KalSeedCollection>     lineH      ; event.getByLabel(line_tag_      , lineH);
    art::Handle<ProtonBunchIntensity>  pbiH       ; event.getByLabel(pbi_tag_       , pbiH);

    TriggerResultsNavigator trigNav(triggerH.product());
    sim_col_        = (simH        .isValid()) ? simH.product()        : nullptr;
    primary_        = (primaryH    .isValid()) ? primaryH.product()    : nullptr;
    mc_digi_col_    = (mc_digiH    .isValid()) ? mc_digiH.product()    : nullptr;
    combo_hits_col_ = (combo_hitsH .isValid()) ? combo_hitsH.product() : nullptr;
    calo_hits_col_  = (calo_hitsH  .isValid()) ? calo_hitsH.product()  : nullptr;
    calo_cluster_mc_col_ = (calo_cluster_mcH.isValid()) ? calo_cluster_mcH.product() : nullptr;
    calo_cluster_mc_assn_ = (calo_cluster_mcassnH.isValid()) ? calo_cluster_mcassnH.product() : nullptr;
    line_col_       = (lineH       .isValid()) ? lineH.product()       : nullptr;
    cluster_col_    = clusterH     .product();
    trig_nav_       = &trigNav;

    if(debug_level_ > 1) std::cout << "[Run1BAna::" << __func__ << "::" << moduleDescription().moduleLabel() << "]"
                                   << " Input from:"
                                   << "\n  " << clusters_tag_.encode().c_str()   << ": N(clusters) = "   << cluster_col_->size()
                                   << "\n  " << calo_hits_tag_.encode().c_str()  << ": N(calo hits) = "  << int((calo_hits_col_)  ? calo_hits_col_->size()  : -1)
                                   << "\n  " << combo_hits_tag_.encode().c_str() << ": N(combo hits) = " << int((combo_hits_col_) ? combo_hits_col_->size() : -1)
                                   << "\n  " << mc_digi_tag_.encode().c_str()    << ": N(MC digis) = "   << int((mc_digi_col_)    ? mc_digi_col_->size()    : -1)
                                   << "\n  " << calo_cluster_mc_tag_.encode().c_str() << ": N(calo cluster MC) = " << int((calo_cluster_mc_col_) ? calo_cluster_mc_col_->size() : -1)
                                   << "\n  " << line_tag_.encode().c_str()       << ": N(lines) = "      << int((line_col_)       ? line_col_->size()       : -1)
                                   << std::endl;

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
    evt_par_.init((pbiH.isValid()) ? pbiH->intensity() : 0);
    const SimParticle* primary_sim = (primary_ && !primary_->primarySimParticles().empty()) ? &(*primary_->primarySimParticles().front()) : nullptr;
    sim_par_.init((primary_sim) ? primary_sim : nullptr, (primary_sim) ? sim_info_[primary_sim->id().asUint()].nhits_ : 0, 0.);

    // If the first event, initialize the trigger path histogram bins for stability
    if(nevt_ == 1) {
      for(int iset = 0; iset < kMaxHists; ++iset) {
        if(!hist_[iset]) continue;
        TH1* h = hist_[iset]->event.trig_paths;
        if(!h) continue;
        for (size_t index = 0; index< trigNav.getTrigPaths().size(); ++index) {
          const std::string path = trigNav.getTrigPathName(index);
          h->GetXaxis()->SetBinLabel(h->FindBin(index), path.c_str());
        }
      }
    }

    //--------------------------------------------------------------------------------------
    // Fill histograms
    //--------------------------------------------------------------------------------------

    const CaloCluster* max_cluster = nullptr;
    const CaloCluster* best_cluster = nullptr;
    for(const auto& cluster : *(cluster_col_)) {
      initClusterPar(cluster_par_, &cluster);



      // Find the highest energy cluster
      if(!max_cluster || max_cluster->energyDep() < cluster.energyDep()) max_cluster = &cluster;

      // Find the "best" cluster in the event
      bool cluster_id = (cluster.energyDep() > 70.
                         && cluster.time() > 600.
                         && cluster.time() < 1650.);
      if(cluster_id) {
        ++evt_par_.n_good_clusters;
        if(!best_cluster || best_cluster->energyDep() < cluster.energyDep()) best_cluster = &cluster;
      }

      // All clusters
      fillHistograms(hist_[1]);

      // Clusters above 70 MeV
      if(cluster.energyDep() > 70.) {
        fillHistograms(hist_[2]);
        if(cluster.diskID() == 0) fillHistograms(hist_[3]);
        else                      fillHistograms(hist_[4]);
      }

      // All clusters passing the cluster ID
      if(cluster_id) fillHistograms(hist_[5]);
    }

    // Per-line (KalSeed) histograms
    if(line_col_) {
      for(const auto& line : *(line_col_)) {
        line_par_.init(&line);
        if(hist_[5]) fillLineHistograms(&hist_[5]->line);
      }
    }

    // All events, highest energy cluster
    initClusterPar(cluster_par_, max_cluster);
    fillHistograms(hist_[0]);

    // above 70 MeV
    if(cluster_par_.cluster && cluster_par_.cluster->energyDep() > 70.) {
      fillHistograms(hist_[10]);
      const float dt = cluster_par_.line_dt();
      if(std::fabs(dt) < 50.) fillHistograms(hist_[11]);
      else                    fillHistograms(hist_[12]);
    }

    // "best" cluster in the event, passing the cluster ID
    if(best_cluster) {
      initClusterPar(cluster_par_, best_cluster);
      fillHistograms(hist_[20]);
      const float dt = cluster_par_.line_dt();
      if(std::fabs(dt) < 50.) fillHistograms(hist_[21]);
      else                    fillHistograms(hist_[22]);
    }


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
