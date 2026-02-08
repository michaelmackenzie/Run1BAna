// CaloEffAna: simple analyzer for primary particle and calorimeter energy
// Created 2026

// framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

// Offline
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"

// ROOT
#include "TH1.h"
#include "TH2.h"

// c++
#include <string>
#include <vector>
#include <iostream>

namespace mu2e {

  class CaloEffAna : public art::EDAnalyzer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> primaryTag     { Name("primary")                 , Comment("PrimaryParticle tag") };
      fhicl::Atom<art::InputTag> caloClusterTag { Name("CaloClusterCollection")   , Comment("CaloCluster collection") };
      fhicl::Atom<art::InputTag> caloShowerTag  { Name("CaloShowerStepCollection"), Comment("CaloShowerStep collection") };
      fhicl::Atom<int>           debugLevel     { Name("debugLevel")              , Comment("Debug level"), 0 };
    };

    // Histograms
    struct Hist_t {
      TH1* h_primary_energy_;
      TH1* h_primary_pdg_;
      TH1* h_nclusters_;
      TH1* h_cluster_energy_;
      TH1* h_max_cluster_energy_;
      TH1* h_total_calo_energy_;
      TH1* h_step_energy_; // individual CaloShowerStep energies
      TH2* h_primary_vs_edep_;
    };
    enum {kMaxHists = 100};

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit CaloEffAna(const Parameters& conf);
    virtual void analyze(const art::Event& event) override;

    void bookHistograms(const int index, const char* title);
    void fillHistograms(Hist_t* Hist);

  private:
    art::InputTag primary_tag_;
    art::InputTag cluster_tag_;
    art::InputTag calo_shower_tag_;
    int debug_level_;

    Hist_t* hists_[kMaxHists];
    const PrimaryParticle*                 primary_     = nullptr;
    const CaloClusterCollection*           cluster_col_ = nullptr;
    const CaloShowerStepCollection*        shower_col_  = nullptr;
  };


  CaloEffAna::CaloEffAna(const Parameters& conf)
    : art::EDAnalyzer{conf}
    , primary_tag_(conf().primaryTag())
    , cluster_tag_(conf().caloClusterTag())
    , calo_shower_tag_(conf().caloShowerTag())
    , debug_level_(conf().debugLevel())
  {
    // register products consumed
    consumes<PrimaryParticle>(primary_tag_);
    consumes<CaloClusterCollection>(cluster_tag_);
    consumes<CaloShowerStepCollection>(calo_shower_tag_);

    for(int i = 0; i < kMaxHists; ++i) hists_[i] = nullptr;
    bookHistograms(0, "all events");
    bookHistograms(1, "edep 1 MeV");
    bookHistograms(2, "edep 10 MeV");
    bookHistograms(3, "edep 50 MeV");
  }

  void CaloEffAna::bookHistograms(const int index, const char* title) {
    if(index >= kMaxHists) throw std::runtime_error("Too many histograms!");
    hists_[index] = new Hist_t;
    auto Hist = hists_[index];

    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory dir = tfs->mkdir(std::format("hist_{}", index), title);

   Hist->h_primary_energy_     = dir.make<TH1F>("primary_energy"    , "Primary energy;Energy (MeV)"              , 150,    0., 150.);
   Hist->h_primary_pdg_        = dir.make<TH1D>("primary_pdg"       , "Primary PDG ID;PDG ID"                    ,  50,  -25.,  25.);
   Hist->h_nclusters_          = dir.make<TH1I>("nclusters"         , "Number of Calo clusters;N clusters"       , 100,    0 , 100 );
   Hist->h_cluster_energy_     = dir.make<TH1F>("cluster_energy"    , "Calo cluster energy;Energy (MeV)"         , 150,    0., 150.);
   Hist->h_max_cluster_energy_ = dir.make<TH1F>("max_cluster_energy", "Max cluster energy;Energy (MeV)"          , 150,    0., 150.);
   Hist->h_total_calo_energy_  = dir.make<TH1F>("total_calo_energy" , "Total Calo energy from steps;Energy (MeV)", 200,    0., 200.);
   Hist->h_step_energy_        = dir.make<TH1F>("calo_step_energy"  , "CaloShowerStep energy;E_{step} (MeV)"     , 100,    0., 100.);
  }

  void CaloEffAna::fillHistograms(Hist_t* Hist) {
    if(!Hist) throw std::runtime_error("Uninitialized histogram book!");

    const SimParticle* primsim = nullptr;
    if(primary_ && !primary_->primarySimParticles().empty()) primsim = &(*primary_->primarySimParticles().front());

    if(primsim) {
      Hist->h_primary_energy_->Fill(primsim->startMomentum().e());
      Hist->h_primary_pdg_->Fill(primsim->pdgId());
      if(debug_level_ > 1) std::cout << "CaloEffAna: primary energy " << primsim->startMomentum().e() << " PDG " << primsim->pdgId() << std::endl;
    }

    if(cluster_col_) {
      Hist->h_nclusters_->Fill(cluster_col_->size());
      const CaloCluster* max_cl = nullptr;
      for(const auto& cl : *cluster_col_) {
        if(!max_cl || cl.energyDep() > max_cl->energyDep()) max_cl = &cl;
        Hist->h_cluster_energy_->Fill(cl.energyDep());
      }
      if(max_cl) Hist->h_max_cluster_energy_->Fill(max_cl->energyDep());
    }

    if(shower_col_) {
      double totalCaloE = 0.;
      for(const auto& css : *shower_col_) {
        const double e = css.energyDepBirks();
        totalCaloE += e;
        Hist->h_step_energy_->Fill(e);
      }
      Hist->h_total_calo_energy_->Fill(totalCaloE);
    }
  }


  void CaloEffAna::analyze(const art::Event& event) {
    // Primary particle
    art::Handle<PrimaryParticle> primaryH;
    event.getByLabel(primary_tag_, primaryH);
    primary_ = (primaryH.isValid()) ? primaryH.product() : nullptr;

    // Calo clusters
    art::Handle<CaloClusterCollection> clusterH;
    event.getByLabel(cluster_tag_, clusterH);
    cluster_col_ = (clusterH.isValid()) ? clusterH.product() : nullptr;

    // Calo shower steps - sum energyDepBirks across the collection
    art::Handle<CaloShowerStepCollection> cssH;
    event.getByLabel(calo_shower_tag_, cssH);
    shower_col_ = (cssH.isValid()) ? cssH.product() : nullptr;

    double totalCaloE = 0.;
    if(shower_col_) {
      for(const auto& css : *shower_col_) {
        const double e = css.energyDepBirks();
        totalCaloE += e;
      }
    }

    // Fill histograms
    fillHistograms(hists_[0]);
    if(totalCaloE >  1.) fillHistograms(hists_[1]);
    if(totalCaloE > 10.) fillHistograms(hists_[2]);
    if(totalCaloE > 50.) fillHistograms(hists_[3]);
  }

}

using mu2e::CaloEffAna;
DEFINE_ART_MODULE(CaloEffAna)
