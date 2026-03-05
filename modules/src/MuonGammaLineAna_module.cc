// ======================================================================
//
// Analyze info for the 1.8 MeV photon line
//
// ======================================================================

// ROOT includes
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "artdaq-core-mu2e/Data/EventHeader.hh"
#include "fhiclcpp/ParameterSet.h"

// Mu2e includes
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"
#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include "Offline/MCDataProducts/inc/CaloShowerSim.hh"
#include "Offline/MCDataProducts/inc/ProcessCode.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/DataProducts/inc/CaloConst.hh"


// std includes
#include <iostream>
#include <memory>
#include <string>
#include <map>

// ======================================================================
namespace mu2e {

  class MuonGammaLineAna : public art::EDAnalyzer {

  public:

  struct CaphriHit_t {
    double energy;
    double time;
    double mc_edep;
    double mc_time;
    int crystal_id;

    CaphriHit_t() { init(); }
    void init(double energy_ = 0., double time_ = 0., double mc_edep_ = 0., double mc_time_ = 0., int crystal_id_ = -1) {
      energy = energy_;
      time = time_;
      mc_edep = mc_edep_;
      mc_time = mc_time_;
      crystal_id = crystal_id_;
    }
  };

  struct Hist_t {
    // observables
    TH1* _hNCaloHits;
    TH1* _hCaloEnergy;
    TH1* _hCaphriEnergy;
    TH1* _hNCaphriHits;
    TH1* _hNTimeClusters;
    TH1* _hNTrackerHits;
    TH1* _hSubRuns;

    // reco info
    TH1* _hNPOT_CaloEnergy;

    // MC info
    TH1* _hGammaEDep;
    TH1* _hHitMCEDep;
    TH1* _hHitRes;
    TH1* _hApproxResp;
    TH1* _hSigFrac;
  };
  enum{kMaxHists = 100};

  struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<std::string> simCol{Name("simCollection"), Comment("SimParticle collection")};
      fhicl::Atom<std::string> caloShowerSimCol{Name("caloShowerSimCollection"), Comment("CaloShowerSim collection")};
      fhicl::Atom<std::string> caloHitMCCol{Name("caloHitMCCollection"), Comment("CaloHitMC collection")};
      fhicl::Atom<std::string> caloHitCol{Name("caloHitCollection"), Comment("CaloHit collection")};
      fhicl::Atom<std::string> caloInfoCol{Name("caloInfo"), Comment("IntensityInfoCalo")};
      fhicl::Atom<std::string> timeClusterInfo{Name("timeClusterInfo"), Comment("IntensityInfoTimeCluster")};
      fhicl::Atom<std::string> trackerHitsInfo{Name("trackerHitsInfo"), Comment("IntensityInfoTrackerHits")};
      fhicl::Atom<std::string> eventHeaderInfo{Name("eventHeaderInfo"), Comment("EventHeader")};
      fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("diagnostic level"), 0};
    };

    explicit MuonGammaLineAna(const art::EDAnalyzer::Table<Config>& config);
    virtual ~MuonGammaLineAna() {}

    virtual void beginJob() override;
    virtual void endJob() override;
    virtual void endSubRun(const art::SubRun& sr) override;

    virtual void analyze(const art::Event& e) override;

    void bookHistograms  (const int index, const std::string name);
    void fillCalo        (mu2e::IntensityInfoCalo const& info       , Hist_t* Hist);
    void fillTrackerHits (mu2e::IntensityInfoTrackerHits const& info, Hist_t* Hist);
    void fillTimeClusters(mu2e::IntensityInfoTimeCluster const& info, Hist_t* Hist);
    void fillEventHeaders(mu2e::EventHeader const& info             , Hist_t* Hist);
    double approximateResponse(double energy);

  private:

    std::string sim_tag_;
    std::string calo_shower_sim_tag_;
    std::string calo_hit_mc_tag_;
    std::string calo_hit_tag_;
    std::string calo_info_tag_;
    std::string time_cluster_info_tag_;
    std::string tracker_hits_info_tag_;
    std::string event_header_info_tag_;
    int _diagLevel;
    TRandom3 _rand;

    size_t _nevents;
    size_t _nsubrun;
    std::map<std::string, size_t> _counter_by_object;
    Hist_t* hists_[kMaxHists];
    int _only_caphri_index = -1;
  };

  // ======================================================================

  MuonGammaLineAna::MuonGammaLineAna(const art::EDAnalyzer::Table<Config>& config) : art::EDAnalyzer{config}
                                                                              , sim_tag_(config().simCol())
                                                                              , calo_shower_sim_tag_(config().caloShowerSimCol())
                                                                              , calo_hit_mc_tag_(config().caloHitMCCol())
                                                                              , calo_hit_tag_(config().caloHitCol())
                                                                              , calo_info_tag_(config().caloInfoCol())
                                                                              , time_cluster_info_tag_(config().timeClusterInfo())
                                                                              , tracker_hits_info_tag_(config().trackerHitsInfo())
                                                                              , event_header_info_tag_(config().eventHeaderInfo())
                                                                              , _diagLevel(config().diagLevel())
                                                                              , _rand(90)
                                                                              , _nevents(0)
                                                                              , _nsubrun(0)
  {
  }

  //--------------------------------------------------------------------------------
  // create the histograms
  //--------------------------------------------------------------------------------
  void MuonGammaLineAna::beginJob() {
    bookHistograms(0, "All events");
    bookHistograms(1, "Signal-like events");
  }

  //--------------------------------------------------------------------------------
  void MuonGammaLineAna::bookHistograms(const int index, const std::string name) {
    if(index >= kMaxHists) throw std::runtime_error("Too many histograms!");
    // book histograms for this index
    hists_[index] = new Hist_t;
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory infoDir = tfs->mkdir(std::format("info_{}", index), name.c_str());

    hists_[index]->_hNCaloHits     = infoDir.make<TH1F>("hNCaloHits", "N(Calorimeter hits);N(hits);Events", 100, 0., 5000.);
    hists_[index]->_hCaloEnergy    = infoDir.make<TH1F>("hCaloEnergy", "Calorimeter energy;Energy (MeV);Events", 100, 0., 6000.);
    hists_[index]->_hCaphriEnergy  = infoDir.make<TH1F>("hCaphriEnergy", "CAPHRI energy;Energy (MeV);Events", 100, 0., 10.);
    hists_[index]->_hNCaphriHits   = infoDir.make<TH1F>("hNCaphriHits", "N(CAPHRI hits);N(hits);Events", 100, 0., 100.);
    hists_[index]->_hNTrackerHits  = infoDir.make<TH1F>("hNTrackerHits", "N(tracker hits);N(hits);Events", 100, 0., 10000.);
    hists_[index]->_hNTimeClusters = infoDir.make<TH1F>("hNTimeClusters", "N(time clusters);N(time clusters);Events", 100, 0., 100.);
    hists_[index]->_hSubRuns       = infoDir.make<TH1F>("hSubRuns", "Sub-Runs;Sub-Run;Events", 100, 0., 100.);
    hists_[index]->_hNPOT_CaloEnergy = infoDir.make<TH1F>("hNPOT_CaloEnergy", "N(POT) from Calorimeter energy;N(POT);Events", 100, 0., 2.e8);
    hists_[index]->_hGammaEDep       = infoDir.make<TH1F>("hGammaEDep", "Gamma energy deposit;Energy (MeV);Events", 100, 0., 10.);
    hists_[index]->_hHitMCEDep       = infoDir.make<TH1F>("hHitMCEDep", "MC energy deposit for hits;Energy (MeV);Events", 100, 0., 10.);
    hists_[index]->_hHitRes          = infoDir.make<TH1F>("hHitRes", "Hit energy residual;Hit energy - MC energy deposited (MeV);Events", 100, -5., 5.);
    hists_[index]->_hApproxResp      = infoDir.make<TH1F>("hApproxResp", "Approximate calorimeter response;Energy (MeV);Events", 100, 0., 10.);
    hists_[index]->_hSigFrac         = infoDir.make<TH1F>("hSigFrac", "Signal energy fraction;Signal energy / Total energy;Events", 110, 0., 1.1);
  }

  void MuonGammaLineAna::endJob() {
    std::cout << "MuonGammaLineAna::" << __func__ << ": Read information for an accumulated " << _nevents << " events\n";
  }

  //--------------------------------------------------------------------------------
  void MuonGammaLineAna::fillCalo(mu2e::IntensityInfoCalo const& info, Hist_t* Hist) {
    Hist->_hNCaloHits  ->Fill(info.nCaloHits());
    Hist->_hCaloEnergy ->Fill(info.caloEnergy());
    Hist->_hNCaphriHits->Fill(info.nCaphriHits());
    for(size_t ihit = 0; ihit < info.nCaphriHits(); ++ihit) {
      if(_only_caphri_index >= 0 && static_cast<int>(ihit) != _only_caphri_index) continue;
      double energy;
      int id;
      info.getCaphriHit(ihit, energy, id);
      Hist->_hCaphriEnergy->Fill(energy);
    }
  }

  //--------------------------------------------------------------------------------
  void MuonGammaLineAna::fillTrackerHits(mu2e::IntensityInfoTrackerHits const& info, Hist_t* Hist) {
     Hist->_hNTrackerHits->Fill(info.nTrackerHits());
    _counter_by_object["tracker"] += 1;
  }

  //--------------------------------------------------------------------------------
  void MuonGammaLineAna::fillTimeClusters(mu2e::IntensityInfoTimeCluster const& info, Hist_t* Hist) {
     Hist->_hNTimeClusters->Fill(info.nProtonTCs());
      _counter_by_object["time_clusters"] += 1;
  }

  //--------------------------------------------------------------------------------
  void MuonGammaLineAna::fillEventHeaders(mu2e::EventHeader const& info, Hist_t* Hist) {
    // Hist->_hSubRuns->Fill(info.subRun());
    _counter_by_object["headers"] += 1;
  }

  //--------------------------------------------------------------------------------
  double MuonGammaLineAna::approximateResponse(double energy) {
    // Approximate calorimeter response function for 1.8 MeV photons, based on a simple Gaussian smearing with a resolution of ~10% at this energy
    const double resolution = 0.1; // 10% resolution
    const double sigma = energy * resolution;
    const double resp = _rand.Gaus(energy, sigma);
    return std::max(0., resp); // truncate at 0
  }

  //--------------------------------------------------------------------------------
  void MuonGammaLineAna::endSubRun(const art::SubRun& sr) {
    _nsubrun += 1;
    hists_[0]->_hSubRuns->Fill(sr.subRun());
  }

  //--------------------------------------------------------------------------------
  void MuonGammaLineAna::analyze(const art::Event& event) {
    // for printout use
    const art::EventNumber_t  eventNumber  = event.event();
    const art::SubRunNumber_t subrunNumber = event.subRun();
    const art::RunNumber_t    runNumber    = event.run();
    _only_caphri_index = -1;

    //---------------------------------------
    // Retrieve the data

    // Calorimeter data
    auto caloHandle = event.getValidHandle<mu2e::IntensityInfoCalo>(calo_info_tag_);
    auto caloInfo = *caloHandle;
    // Tracker hit data
    auto trkHandle = event.getValidHandle<mu2e::IntensityInfoTrackerHits>(tracker_hits_info_tag_);
    auto trkInfo = *trkHandle;
    // Time cluster data
    auto tcHandle = event.getValidHandle<mu2e::IntensityInfoTimeCluster>(time_cluster_info_tag_);
    auto tcInfo = *tcHandle;
    // Event header data
    auto headerHandle = event.getValidHandle<mu2e::EventHeader>(event_header_info_tag_);
    auto headerInfo = *headerHandle;
    // SimParticle data
    auto simHandle = event.getValidHandle<mu2e::SimParticleCollection>(sim_tag_);
    auto simCol = *simHandle;
    // Shower info
    auto caloShowerSimHandle = event.getValidHandle<mu2e::CaloShowerSimCollection>(calo_shower_sim_tag_);
    auto caloShowerSimCol = *caloShowerSimHandle;
    // Calo hits
    auto caloHitHandle = event.getValidHandle<mu2e::CaloHitCollection>(calo_hit_tag_);
    auto caloHitCol = *caloHitHandle;
    // Calo MC hits
    auto caloMCHitsHandle = event.getValidHandle<CaloHitMCCollection>(calo_hit_mc_tag_);
    auto caloMCHitCol = *caloMCHitsHandle;

    if(_diagLevel > 1) {
      std::cout << "MuonGammaLineAna::" << __func__ << ": Subrun "
                << runNumber << ":" << subrunNumber << ":" << eventNumber
                << ": SimParticles = " << simCol.size()
                << std::endl;
    }

    if(caloMCHitCol.size() != caloHitCol.size()) {
      throw std::runtime_error("CaloHit and CaloHitMC collections have different sizes!");
    }
    fillCalo        (caloInfo  , hists_[0]);
    fillTrackerHits (trkInfo   , hists_[0]);
    fillTimeClusters(tcInfo    , hists_[0]);
    fillEventHeaders(headerInfo, hists_[0]);

    // for(const auto& sim_pair : simCol) {
    //   const auto& sim = sim_pair.second;
    //   if(sim.pdgId() == 22 && sim.startMomentum().e() > 1.805 && sim.startMomentum().e() < 1.815 && sim.creationCode() == mu2e::ProcessCode::mu2eMuonCaptureAtRest) {
    //     if(_diagLevel > 0) {
    //       std::cout << "Found sim particle with PDG = " << sim.pdgId() << " and E = " << sim.startMomentum().e() << " MeV" << std::endl;
    //     }
    //     fillCalo        (caloInfo  , hists_[1]);
    //     fillTrackerHits (trkInfo   , hists_[1]);
    //     fillTimeClusters(tcInfo    , hists_[1]);
    //     fillEventHeaders(headerInfo, hists_[1]);
    //     break;
    //   }
    // }

    // loop through the hit collection
    for(size_t hit_index = 0; hit_index < caloHitCol.size(); ++hit_index) {
      const auto& hit = caloHitCol[hit_index];
      if(hit.energyDep() < 1. || hit.energyDep() > 5.) continue; // only look at hits with 1 < energy < 5 MeV
      const auto& mc_hit = caloMCHitCol[hit_index];
      const double mc_edep = mc_hit.totalEnergyDep();
      hists_[0]->_hApproxResp->Fill(approximateResponse(mc_edep));
      hists_[0]->_hHitMCEDep->Fill(mc_edep);
      hists_[0]->_hHitRes->Fill(hit.energyDep() - mc_edep);
      double signal_edep = 0.;
      for(const auto& edep : mc_hit.energyDeposits()) {
        if(edep.sim().isNull()) continue;
        const auto sim = &(*edep.sim());
        const bool is_signal = sim->pdgId() == 22 &&
                               sim->startMomentum().e() > 1.805 &&
                               sim->startMomentum().e() < 1.815 &&
                               sim->creationCode() == mu2e::ProcessCode::mu2eMuonCaptureAtRest;
        if(is_signal) signal_edep += edep.energyDep();
      }
      hists_[0]->_hSigFrac->Fill(signal_edep / mc_edep);

      if(signal_edep > 0.5*mc_edep) { // mostly a signal hit
        hists_[1]->_hApproxResp->Fill(approximateResponse(mc_edep));
        hists_[1]->_hHitMCEDep->Fill(mc_edep);
        hists_[1]->_hHitRes->Fill(hit.energyDep() - mc_edep);
        hists_[1]->_hSigFrac->Fill(signal_edep / mc_edep);
      }
    }

    // Loop through shower edep
    for(const auto& shower : caloShowerSimCol) {
      if(shower.sim().isNull()) continue;
      const auto shower_sim = &(*shower.sim());
      const bool caphri_crystal = std::find(std::begin(mu2e::CaloConst::_caphriId), std::end(mu2e::CaloConst::_caphriId), shower.crystalID()) != std::end(mu2e::CaloConst::_caphriId);

      // if(caphri_crystal) {
      //   const double edep = shower.energyDep();
      //   const double resp = approximateResponse(edep);
      //   hists_[0]->_hApproxResp->Fill(resp);
      // }

      // Check if this is a CAPHRI hit by the signal photon
      if(caphri_crystal && shower_sim->pdgId() == 22
         && shower_sim->startMomentum().e() > 1.805 && shower_sim->startMomentum().e() < 1.815
         && shower_sim->creationCode() == mu2e::ProcessCode::mu2eMuonCaptureAtRest) {

        // Check for the corresponding CaloHit
        int matching_hit_index = -1;
        for(size_t i = 0; i < caloHitCol.size(); ++i) {
          const auto& hit = caloHitCol[i];
          // std::cout << "Shower edep = " << shower.energyDep() << " MeV, crystal ID = " << shower.crystalID() << ", time = " << shower.time() << " ns"
          //           << " CaloHit energy = " << hit.energyDep() << " MeV, time = " << hit.time() << " ns" << std::endl;
          if(hit.crystalID() == shower.crystalID() && std::abs(hit.time() - shower.time()) < 100.) { // 100 ns time window
            matching_hit_index = i;
            break;
          }
        }
        if(matching_hit_index == -1) {
          if(_diagLevel > 0 && shower.energyDep() > 1.) { // only print if large deposit
            std::cout << "No matching CaloHit found for CAPHRI shower with crystal ID = " << shower.crystalID() << " and time = " << shower.time() << std::endl;
          }
          continue;
        }
        const auto& matching_hit = caloHitCol[matching_hit_index];
        const auto& mc_hit = caloMCHitCol[matching_hit_index];
        if(_diagLevel > 0) {
          std::cout << event.id() << " : Found CAPHRI CaloShowerSim with PDG = " << shower_sim->pdgId()
                    << " and E = " << shower_sim->startMomentum().e() << " MeV"
                    << " edep = " << shower.energyDep() << " MeV"
                    << " crystal ID = " << shower.crystalID()
                    << " time = " << shower.time() << " ns"
                    << " and matching CaloHit with energy = " << matching_hit.energyDep() << " MeV and time = " << matching_hit.time() << " ns"
                    << " MC edep = " << mc_hit.totalEnergyDep() << " MeV"
                    << " N(CAPHRI hits) = " << caloInfo.nCaphriHits()
                    << std::endl;

        }
        double energy;
        int id;
        bool crystal_found = false;
        for(size_t ihit = 0; ihit < caloInfo.nCaphriHits(); ++ihit) {
          caloInfo.getCaphriHit(ihit, energy, id);
          if(_diagLevel > 0) std::cout << "  CAPHRI hit " << ihit << ": energy = " << energy << " MeV, ID = " << id << std::endl;
          if(id == shower.crystalID() && std::abs(energy - matching_hit.energyDep()) < 0.2) { // check if this CAPHRI hit matches the CaloHit
            if(crystal_found && _diagLevel > 0) {
              std::cout << "Warning: Multiple CAPHRI hits found matching the same shower and CaloHit!" << std::endl;
            }
            crystal_found = true;
            _only_caphri_index = ihit;
            fillCalo        (caloInfo  , hists_[1]);
            fillTrackerHits (trkInfo   , hists_[1]);
            fillTimeClusters(tcInfo    , hists_[1]);
            fillEventHeaders(headerInfo, hists_[1]);
            hists_[1]->_hGammaEDep->Fill(shower.energyDep());
            // hists_[1]->_hHitRes->Fill(energy - shower.energyDep());
            // hists_[1]->_hApproxResp->Fill(approximateResponse(shower.energyDep()));
          }
        }
        if(crystal_found) {
          break;
        }
      }
    }
    _only_caphri_index = -1;
     _nevents += 1;
  }

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::MuonGammaLineAna)
