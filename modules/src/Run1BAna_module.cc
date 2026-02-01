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
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"

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
  class Run1BAna : public art::EDAnalyzer
  {
  public:

    //--------------------------------------------------------------------------------------
    // Input config
    //--------------------------------------------------------------------------------------

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>      simCol        { Name("simCollection")        , Comment("simCollection")                       };
      fhicl::Atom<art::InputTag>      digiCol       { Name("strawDigiMCCollection"), Comment("strawDigiMCCollection")               };
      fhicl::Atom<art::InputTag>      caloClusterCol{ Name("caloClusterCollection"), Comment("caloClusterCollection")               };
      fhicl::Atom<std::string>        trigProcess   { Name("triggerProcess")       , Comment("Process name for the trigger results")};
      fhicl::Atom<art::InputTag>      pbi           { Name("PBI")                  , Comment("ProtonBunchIntensity tag")            };
      fhicl::Atom<int>                debugLevel    { Name("debugLevel")           , Comment("debugLevel")                         , 0 };
    };


    //--------------------------------------------------------------------------------------
    // Histograms
    //--------------------------------------------------------------------------------------

    // Per event
    struct EventHist_t {
      TH1* npot;
      TH1* n_mc_digis;
      TH1* nhits;
      TH1* nclusters;

      TH1* trig_bits;
      TH1* trig_paths;
    };

    // Per cluster
    struct ClusterHist_t {
      TH1* energy;
      TH1* time;
      TH1* radius;
      TH1* ncr;
      TH1* disk;
    };

    // Total
    struct Hist_t {
      EventHist_t   event;
      ClusterHist_t cluster;
    };

    //--------------------------------------------------------------------------------------
    // Internal data structures
    //--------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------
    struct EventPar_t {
      long long npot;

      EventPar_t() {
        init();
      }

      void init(long long np = 0) {
        npot = np;
      }
    };

    //--------------------------------------------------------------------------------------
    struct ClusterPar_t {
      const CaloCluster* cluster;
      float r;

      ClusterPar_t() {
        init();
      }

      void init(const CaloCluster* cl = nullptr) {
        cluster = cl;
        r = 0.f;
        if(!cl) return;

        const float x = cluster->cog3Vector().x();
        const float y = cluster->cog3Vector().y();
        r = std::sqrt(x*x + y*y);
      }
    };

    //--------------------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------------------

    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit Run1BAna(const Parameters& conf);
    virtual void analyze(const art::Event& event) override;
    virtual void beginJob();
    virtual void beginRun(const art::Run&   run   );
    virtual void endRun(const art::Run& run ) override;

  private:

    void bookClusterHistograms(ClusterHist_t* Hist, const char* name);
    void bookEventHistograms(EventHist_t* Hist, const char* name);
    void bookHistograms(const int index, const char* name);
    void fillClusterHistograms(ClusterHist_t* Hist);
    void fillEventHistograms(EventHist_t* Hist);
    void fillHistograms(Hist_t* Hist);


    //--------------------------------------------------------------------------------------
    // Data
    //--------------------------------------------------------------------------------------

    art::InputTag  _simTag;
    art::InputTag  _digTag;
    art::InputTag  _clsTag;
    art::InputTag  _trgTag;
    art::InputTag  _pbiTag;
    int            _debugLevel;

    enum {kMaxHists = 100};
    Hist_t* _hist[kMaxHists];
    const SimParticleCollection*    _simcol = nullptr;
    const StrawDigiMCCollection*    _digcol = nullptr;
    const CaloClusterCollection*    _clscol = nullptr;
    const TriggerResultsNavigator*  _trgnav = nullptr;
    EventPar_t                      _evt_par;
    ClusterPar_t                    _cluster_par;

    unsigned  long                  _nevt;
  };

  //--------------------------------------------------------------------------------------
  Run1BAna::Run1BAna(const Parameters& config):
    art::EDAnalyzer{config}
    , _simTag            (config().simCol())
    , _digTag            (config().digiCol())
    , _clsTag            (config().caloClusterCol())
    , _trgTag            ("TriggerResults::" + config().trigProcess())
    , _pbiTag            (config().pbi())
    , _debugLevel        (config().debugLevel())
    , _nevt              (0)
    {
      for(int i = 0; i < kMaxHists; ++i) _hist[i] = nullptr;
      bookHistograms(0, "all_events");
      bookHistograms(1, "all_clusters");
      bookHistograms(2, "70MeV_clusters");
      bookHistograms(3, "70MeV_disk0");
      bookHistograms(4, "70MeV_disk1");

      bookHistograms(10, "max_70MeV");
    }

  //--------------------------------------------------------------------------------------
  void Run1BAna::beginJob() {
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::beginRun(const art::Run & run){
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookEventHistograms(EventHist_t* Hist, const char* name) {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(name);
    Hist->npot          = dir.make<TH1F>("npot"         , "N(POT);"                     ,  100,    0.,   1.e7);
    Hist->n_mc_digis    = dir.make<TH1F>("n_mc_digis"   , "N(MC digis);"                ,   25,    0., 10000.);
    Hist->nhits         = dir.make<TH1F>("nhits"        , "N(tracker hits);"            ,   25,    0., 10000.);
    Hist->nclusters     = dir.make<TH1D>("nclusters"    , "N(calo clusters);"           ,  100,    0.,   100.);

    Hist->trig_bits     = dir.make<TH1D>("trig_bits"    , "Trigger bits;"               , 1000,    0.,   1000);
    Hist->trig_paths    = dir.make<TH1D>("trig_paths"   , "Trigger paths;"              ,  100,    0.,    100);

  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookClusterHistograms(ClusterHist_t* Hist, const char* name) {
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(name);
    Hist->energy        = dir.make<TH1F>("energy"  , "Cluster energy;Energy (MeV);"     ,  100,    0.,    300.);
    Hist->time          = dir.make<TH1F>("time"    , "Cluster time;Time (ns);"          ,  100,    0.,   2000.);
    Hist->radius        = dir.make<TH1F>("radius"  , "Cluster radial position;R (mm);"  ,  100,    0.,    700.);
    Hist->ncr           = dir.make<TH1D>("ncr"     , "Number of crystals;N(Crystals);"  ,   20,    0.,     20.);
    Hist->disk          = dir.make<TH1D>("disk"    , "Disk ID;Disk ID;"                 ,    3,   -1.,      2.);
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::bookHistograms(const int index, const char* name) {
    if(index >= kMaxHists) throw std::runtime_error("Too many histograms!");
    _hist[index] = new Hist_t;
    auto Hist = _hist[index];
    bookEventHistograms  (&Hist->event  , std::format("evt_{:s}", name).c_str());
    bookClusterHistograms(&Hist->cluster, std::format("cls_{:s}", name).c_str());
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillClusterHistograms(ClusterHist_t* Hist) {
    if(!Hist) return;
    auto Cluster = _cluster_par.cluster;
    if(!Cluster) return;

    Hist->energy->Fill(Cluster->energyDep());
    Hist->time->Fill(Cluster->time());
    Hist->radius->Fill(_cluster_par.r);
    Hist->ncr->Fill(Cluster->caloHitsPtrVector().size());
    Hist->disk->Fill(Cluster->diskID());
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillEventHistograms(EventHist_t* Hist) {
    if(!Hist) return;
    Hist->npot          ->Fill(_evt_par.npot);
    Hist->n_mc_digis    ->Fill((_digcol) ? _digcol->size() : 0);
    Hist->nhits         ->Fill(0);
    Hist->nclusters     ->Fill((_clscol) ? _clscol->size() : 0);

    // Trigger information
    for (size_t index = 0; index < _trgnav->getTrigPaths().size(); ++index) {
      const std::string path = _trgnav->getTrigPathName(index);
      if(_trgnav->accepted(path)) {
        Hist->trig_bits ->Fill(_trgnav->findTrigPathID(path));
        Hist->trig_paths->Fill(path.c_str(), 1.);
      }
    }
  }

  //--------------------------------------------------------------------------------------
  void Run1BAna::fillHistograms(Hist_t* Hist) {
    if(!Hist) return;
    fillEventHistograms(&Hist->event);
    fillClusterHistograms(&Hist->cluster);
  }


  //--------------------------------------------------------------------------------------
  void Run1BAna::analyze(const art::Event& event){
    ++_nevt;

    //--------------------------------------------------------------------------------------
    // Retrieve the collections
    //--------------------------------------------------------------------------------------

    // auto simH = event.getValidHandle<SimParticleCollection>(_simTag);
    // auto digH = event.getValidHandle<StrawDigiMCCollection>(_digTag);
    auto clsH = event.getValidHandle<CaloClusterCollection>(_clsTag);
    auto trgH = event.getValidHandle<art::TriggerResults>  (_trgTag);
    auto pbiH = event.getValidHandle<ProtonBunchIntensity> (_pbiTag);

    TriggerResultsNavigator trigNav(trgH.product());
    // _simcol = simH.product();
    // _digcol = digH.product();
    _clscol = clsH.product();
    _trgnav = &trigNav;
    // art::Handle<mu2e::ProtonBunchIntensity> pbiHandle;
    // AnEvent->getByLabel("PBISim", pbiHandle);


    if(_debugLevel > 1) std::cout << "[Run1BAna::" << __func__ << "::" << moduleDescription().moduleLabel() << "]"
                                  << " Input from:"
                                  // << " " << _simTag.encode().c_str() << " N(sims) = "<< _simcol->size()
                                  // << " " << _digTag.encode().c_str() << " N(digis) = "<< _digcol->size()
                                  << " " << _clsTag.encode().c_str() << " N(clusters) = " << _clscol->size()
                                  << std::endl;

    //--------------------------------------------------------------------------------------
    // Initialize fields
    //--------------------------------------------------------------------------------------

    _cluster_par.init();
    _evt_par.init(pbiH->intensity());

    // If the first event, initialize the trigger path histogram bins for stability
    if(_nevt == 1) {
      for(int iset = 0; iset < kMaxHists; ++iset) {
        if(!_hist[iset]) continue;
        TH1* h = _hist[iset]->event.trig_bits;
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
    for(const auto& cluster : *(_clscol)) {
      _cluster_par.init(&cluster);

      // Find the highest energy cluster
      if(!max_cluster || max_cluster->energyDep() < cluster.energyDep()) max_cluster = &cluster;

      // All clusters
      fillHistograms(_hist[1]);

      // Clusters above 70 MeV
      if(cluster.energyDep() > 70.) {
        fillHistograms(_hist[2]);
        if(cluster.diskID() == 0) fillHistograms(_hist[3]);
        else                      fillHistograms(_hist[4]);
      }
    }

    // All events, highest energy cluster
    _cluster_par.init(max_cluster);
    fillHistograms(_hist[0]);
    if(_cluster_par.cluster && _cluster_par.cluster->energyDep() > 70.) {
      fillHistograms(_hist[10]);
    }
  }


  //--------------------------------------------------------------------------------------
  void Run1BAna::endRun(const art::Run& run) {
  }

}
using mu2e::Run1BAna;
DEFINE_ART_MODULE(Run1BAna)
