// ======================================================================
//
// Get RPC normalization info
//
// ======================================================================

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
#include "Offline/MCDataProducts/inc/GenEventCount.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/Mu2eUtilities/inc/SimParticleGetTau.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

// C++ includes
#include <iostream>
#include <memory>
#include <string>

// ROOT includes
#include "TH1.h"
#include "TH2.h"

// ======================================================================
namespace mu2e {

  class RPCCounting : public art::EDAnalyzer {

  public:

  struct Hist_t {
    TH1* pionTime;
    TH1* pionMomentum;
    TH1* pionZ;
    TH2* pionXY;
  };
  enum{kMaxHists = 100};

  struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<std::string> simCol{Name("simCollection"), Comment("SimParticle collection")};
      fhicl::Atom<std::string> genEventCount{Name("genEventCount"), Comment("GenEventCount")};
      fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("diagnostic level"), 0};
    };

    explicit RPCCounting(const art::EDAnalyzer::Table<Config>& config);
    virtual ~RPCCounting() {}

    virtual void beginJob() override;
    virtual void endJob() override;
    virtual void beginSubRun(const art::SubRun& sr) override;

    virtual void analyze(const art::Event& e) override;

    void bookHistograms  (const int index, const std::string name);
    void fillHistograms  (const int index, const SimParticle& sim, const double weight = 1.);

  private:

    std::string sim_tag_;
    std::string gen_event_count_tag_;
    int diag_level_;

    size_t nevents_;
    size_t ngen_;
    double sum_wt_;
    Hist_t* hists_[kMaxHists];
  };

  // ======================================================================

  RPCCounting::RPCCounting(const art::EDAnalyzer::Table<Config>& config) : art::EDAnalyzer{config}
                                                                         , sim_tag_(config().simCol())
                                                                         , gen_event_count_tag_(config().genEventCount())
                                                                         , diag_level_(config().diagLevel())
                                                                         , nevents_(0)
                                                                         , ngen_(0)
                                                                         , sum_wt_(0.)
  {
  }

  //--------------------------------------------------------------------------------
  // create the histograms
  //--------------------------------------------------------------------------------
  void RPCCounting::beginJob() {
    for(int i=0; i<kMaxHists; ++i) hists_[i] = nullptr;
    bookHistograms(0, "No weight");
    bookHistograms(1, "RPC weight");
  }

  //--------------------------------------------------------------------------------
  void RPCCounting::endJob() {
    std::cout << "RPCCounting::" << __func__ << ": Read information for an accumulated " << nevents_ << " events\n";
    std::cout << "RPCCounting::" << __func__ << ": Total GenEventCount = " << ngen_ << std::endl;
    std::cout << "RPCCounting::" << __func__ << ": Total weight = " << sum_wt_ << std::endl;
  }

  //--------------------------------------------------------------------------------
  void RPCCounting::beginSubRun(const art::SubRun& sr) {
    // get the GenEventCount object for this subrun
    auto genEventCountHandle = sr.getValidHandle<mu2e::GenEventCount>(gen_event_count_tag_);
    ngen_ += genEventCountHandle->count();
    if(diag_level_ > 0) {
      std::cout << "RPCCounting::" << __func__ << ": Subrun " << sr.run() << ":" << sr.subRun()
                << ": Total GenEventCount = " << ngen_ << std::endl;
    }
  }

  //--------------------------------------------------------------------------------
  void RPCCounting::bookHistograms(const int index, const std::string name) {
    if(index >= kMaxHists) throw std::runtime_error("Too many histograms!");
    // book histograms for this index
    hists_[index] = new Hist_t;
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory infoDir = tfs->mkdir(std::format("info_{}", index), name.c_str());

    hists_[index]->pionTime        = infoDir.make<TH1F>("pionTime", "Pion stop time;time (ns);", 500, 0., 2000.);
    hists_[index]->pionZ           = infoDir.make<TH1F>("pionZ", "Pion stop z;z (cm);", 500, 0., 10000.);
    hists_[index]->pionXY          = infoDir.make<TH2F>("pionXY", "Pion stop x vs y;x (mm);y (mm);", 500, -1000., 1000., 500, -1000., 1000.);
  }

  //--------------------------------------------------------------------------------
  void RPCCounting::fillHistograms(const int index, const SimParticle& sim, const double weight) {
    if(index >= kMaxHists) throw std::runtime_error("Too many histograms!");
    if(!hists_[index]) throw std::runtime_error("Histogram not booked!");
    hists_[index]->pionTime->Fill(sim.endGlobalTime(), weight);
    hists_[index]->pionZ->Fill(sim.endPosition().z(), weight);
    hists_[index]->pionXY->Fill(sim.endPosition().x()+3904., sim.endPosition().y(), weight);
  }

  //--------------------------------------------------------------------------------
  void RPCCounting::analyze(const art::Event& event) {
    ++nevents_;

    // for printout use
    const art::EventNumber_t  eventNumber  = event.event();
    const art::SubRunNumber_t subrunNumber = event.subRun();
    const art::RunNumber_t    runNumber    = event.run();

    //---------------------------------------
    // Retrieve the data

    // sim particles
    auto simHandle = event.getValidHandle<mu2e::SimParticleCollection>(sim_tag_);
    auto simCol = *simHandle;

    if(diag_level_ > 1) {
      std::cout << "RPCCounting::" << __func__ << ": Subrun "
                << runNumber << ":" << subrunNumber << ":" << eventNumber
                << ": SimParticles = " << simCol.size()
                << std::endl;
    }

    // loop through the pion collection and fill histograms
    for(const auto& [id, sim] : simCol) {
      if(std::abs(sim.pdgId()) == 211 && sim.stoppingCode() == ProcessCode::mu2eKillerVolume) {
        const auto& gc = *GlobalConstantsHandle<PhysicsParams>();
        const double tau = SimParticleGetTau::calculate(sim, {PDGCode::pi_plus, PDGCode::pi_minus}, gc);
        const double proper_time = sim.endProperTime();
        const double weight = std::exp(-tau);
        fillHistograms(0, sim);
        fillHistograms(1, sim, weight);
        sum_wt_ += weight;
        if(diag_level_ > 2) {
          std::cout << "RPCCounting::" << __func__ << event.id()
                    << ": Found pion stop: time = " << sim.endGlobalTime()
                    << ", process code =" << sim.stoppingCode()
                    << ", proper time = " << proper_time
                    << ", tau = " << tau
                    << ", weight = " << weight
                    << std::endl;
        }
      }
    }
  }

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::RPCCounting)
