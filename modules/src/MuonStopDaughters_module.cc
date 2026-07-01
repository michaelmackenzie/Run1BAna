//
//  Measure G4 information for muon daughters
//  Michael MacKenzie, 2026
//
// framework
#include "fhiclcpp/types/Atom.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/TriggerResults.h"

// Mu2e
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TString.h"

// c++
#include <array>
#include <cmath>
#include <set>
#include <string>
#include <vector>
#include <iostream>

// local
#include "Run1BAna/modules/inc/SimUtils.hh"

using namespace CLHEP;

namespace mu2e
{
  class MuonStopDaughters : public art::EDAnalyzer
  {
  public:

    //--------------------------------------------------------------------------------------
    // Input config
    //--------------------------------------------------------------------------------------

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>      simCol           { Name("simCollection")          , Comment("Sim partcile collection") };
      fhicl::Atom<art::InputTag>      physVolInfoInput { Name("physVolInfoInput")       , Comment("PhysicalVolumeInfoMultiCollection input") };
      fhicl::Atom<bool>               useEventLevelVolumeInfo { Name("useEventLevelVolumeInfo"), Comment("Get PhysicalVolumeInfoMultiCollection from Event instead of SubRun"), false };
      fhicl::Atom<std::string>        material         { Name("material")               , Comment("Stop material")           };
      fhicl::Atom<int>                debugLevel       { Name("debugLevel")             , Comment("debugLevel")              , 0 };
    };


    // Histograms
    struct Hist_t {
      TH1* start_energy_ = nullptr;
      TH2* start_xy_ = nullptr;
      TH2* start_rz_ = nullptr;
    };

    //--------------------------------------------------------------------------------------
    // Methods
    //--------------------------------------------------------------------------------------

    using Parameters = art::EDAnalyzer::Table<Config>;

    explicit MuonStopDaughters(const Parameters& conf);
    virtual void analyze(const art::Event& event) override;
    virtual void beginJob();
    virtual void endJob();
    virtual void beginRun(const art::Run&   run   );
    virtual void beginSubRun(const art::SubRun& subrun);
    virtual void endRun(const art::Run& run ) override;

  private:

    void bookHistograms(const int index, const char* name);
    void fillHistograms(Hist_t* Hist, const SimParticle& particle);
    bool isTrackedParticle(PDGCode::type pdgId, size_t& index) const;
    bool acceptedProcess(ProcessCode code) const;

    template <class PRINCIPAL> void initVols(const PRINCIPAL& principal);
    //--------------------------------------------------------------------------------------
    // Data
    //--------------------------------------------------------------------------------------

    art::InputTag  sim_tag_;
    art::InputTag  phys_vol_info_tag_;
    bool           use_event_level_volume_info_;
    std::string    material_;
    int            debug_level_;

    enum {kMaxHists = 100};
    Hist_t* hist_[kMaxHists];
    TH1* hist_norm_;
    const PhysicalVolumeInfoMultiCollection* vols_ = nullptr;
    const SimParticleCollection*           sim_col_ = nullptr;
    const art::Event*                      event_ = nullptr;
    unsigned  long                         nevt_, ngen_;
    unsigned long                          n_muon_stops_;
    std::array<unsigned long, 5>           daughter_counts_;
  };

  //--------------------------------------------------------------------------------------
  MuonStopDaughters::MuonStopDaughters(const Parameters& config):
    art::EDAnalyzer{config}
    , sim_tag_            (config().simCol())
    , phys_vol_info_tag_  (config().physVolInfoInput())
    , use_event_level_volume_info_(config().useEventLevelVolumeInfo())
    , material_           (config().material())
    , debug_level_        (config().debugLevel())
    , nevt_               (0)
    , ngen_               (0)
    , n_muon_stops_       (0)
    , daughter_counts_    {0, 0, 0, 0, 0}
  {
    // Normalization histogram
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir("data", "Data");
    hist_norm_ = dir.make<TH1D>("norm", "Normalization counts", 1, 0., 1.);

    // Analysis histograms
    for(int i = 0; i < kMaxHists; ++i) hist_[i] = nullptr;
    bookHistograms(0, "all_events");
    bookHistograms(1, "protons");
    bookHistograms(2, "neutrons");
    bookHistograms(3, "electrons");
    bookHistograms(4, "photons");
    bookHistograms(5, "deuteron");
  }

  //--------------------------------------------------------------------------------------
  void MuonStopDaughters::beginJob() {
  }

  //--------------------------------------------------------------------------------------
  void MuonStopDaughters::beginRun(const art::Run & run) {
  }

  //--------------------------------------------------------------------------------------
  void MuonStopDaughters::beginSubRun(const art::SubRun& subrun) {
    if(!use_event_level_volume_info_) initVols(subrun);
  }

  //--------------------------------------------------------------------------------------
  void MuonStopDaughters::bookHistograms(const int index, const char* title) {
    if(index >= kMaxHists) throw std::runtime_error("Too many histograms!");
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(title, title);
    hist_[index] = new Hist_t;
    hist_[index]->start_energy_ = dir.make<TH1D>("start_energy", "Daughter start kinetic energy;Kinetic energy [MeV];Entries", 200, 0., 200.);
    hist_[index]->start_xy_ = dir.make<TH2D>("start_xy", "Daughter start position;x [mm];y [mm]", 200, -5000., 5000., 200, -5000., 5000.);
    hist_[index]->start_rz_ = dir.make<TH2D>("start_rz", "Daughter start position;r [mm];z [mm]", 200, 0., 5000., 300, -15000., 15000.);
  }

  //--------------------------------------------------------------------------------------
  void MuonStopDaughters::fillHistograms(Hist_t* Hist, const SimParticle& particle) {
    if(!Hist) return;

    const auto position = particle.startPosition();
    const auto momentum = particle.startMomentum();
    const double mass = momentum.m();
    const double kinetic_energy = momentum.e() - mass;
    const double radius = std::hypot(position.x(), position.y());

    Hist->start_energy_->Fill(kinetic_energy);
    Hist->start_xy_->Fill(position.x(), position.y());
    Hist->start_rz_->Fill(radius, position.z());
  }

  //--------------------------------------------------------------------------------------
  bool MuonStopDaughters::isTrackedParticle(PDGCode::type pdgId, size_t& index) const {
    switch(pdgId) {
      case PDGCode::proton:  index = 0; return true;
      case PDGCode::neutron: index = 1; return true;
      case PDGCode::e_minus:
      case PDGCode::e_plus:  index = 2; return true;
      case PDGCode::gamma:   index = 3; return true;
      case PDGCode::deuteron:index = 4; return true;
      default: return false;
    }
  }

  //--------------------------------------------------------------------------------------
  bool MuonStopDaughters::acceptedProcess(ProcessCode code) const {
    switch(code) {
    case ProcessCode::muIoni: return false;
    case ProcessCode::EMCascade: return false;
    default: return true;
    }
  }

  //--------------------------------------------------------------------------------------
  template <class PRINCIPAL>
  void MuonStopDaughters::initVols(const PRINCIPAL& principal) {
    const auto& volume_handle = principal.template getValidHandle<PhysicalVolumeInfoMultiCollection>(phys_vol_info_tag_);
    vols_ = &*volume_handle;
  }

  //--------------------------------------------------------------------------------------
  void MuonStopDaughters::analyze(const art::Event& event){
    ++nevt_;
    hist_norm_->Fill(0.);
    event_ = &event;

    if(use_event_level_volume_info_) initVols(event);

    //--------------------------------------------------------------------------------------
    // Retrieve the collections
    //--------------------------------------------------------------------------------------

    art::Handle<SimParticleCollection> simH                  ; event.getByLabel(sim_tag_            , simH);
    sim_col_              = (simH                .isValid()) ? simH.product()                 : nullptr;

    //--------------------------------------------------------------------------------------
    // Analyze the data
    //--------------------------------------------------------------------------------------

    if(!sim_col_ || !vols_) return;

    PhysicalVolumeMultiHelper vi(vols_);
    std::set<SimParticle::key_type> counted_muons;
    for(const auto& sim : *sim_col_) {
      const SimParticle& particle = sim.second;
      if(!particle.parent().isNonnull() || !particle.parent().isAvailable()) continue;

      const SimParticle& parent = *particle.parent();
      if(parent.pdgId() != PDGCode::mu_minus) continue;

      const auto& parent_stop_material = vi.endVolume(parent).materialName();
      const auto& daughter_start_material = vi.startVolume(particle).materialName();
      if(parent_stop_material != material_ || daughter_start_material != material_) continue;
      if(!acceptedProcess(particle.creationCode())) continue;
      if(debug_level_ > 0) {
        std::cout << "[MuonStopDaughter::" << __func__ << "] Found daughter: pdg = "
                  << particle.pdgId() << " process = "
                  << particle.creationCode()
                  << std::endl;
      }

      counted_muons.insert(parent.id());

      fillHistograms(hist_[0], particle);

      size_t particle_index = 0;
      if(isTrackedParticle(particle.pdgId(), particle_index)) {
        ++daughter_counts_[particle_index];
        fillHistograms(hist_[particle_index + 1], particle);
      }

    }

    n_muon_stops_ += counted_muons.size();

  } // end analyze


  //--------------------------------------------------------------------------------------
  void MuonStopDaughters::endRun(const art::Run& run) {
  }

  //--------------------------------------------------------------------------------------
  void MuonStopDaughters::endJob() {
    static const std::array<const char*, 5> labels = {{"proton", "neutron", "electron", "photon", "deuteron"}};
    std::cout
      << "MuonStopDaughters stats: muon stops in " << material_ << " = " << n_muon_stops_
      << std::endl;

    for(size_t index = 0; index < labels.size(); ++index) {
      const double rate = (n_muon_stops_ > 0) ? static_cast<double>(daughter_counts_[index]) / n_muon_stops_ : 0.;
      std::cout
        << "MuonStopDaughters production rate: " << labels[index]
        << " count = " << daughter_counts_[index]
        << ", per muon stop = " << rate
        << std::endl;
    }
  }

}
using mu2e::MuonStopDaughters;
DEFINE_ART_MODULE(MuonStopDaughters)
