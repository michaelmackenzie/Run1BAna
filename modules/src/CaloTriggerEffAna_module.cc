#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/TriggerResults.h"
#include "art_root_io/TFileService.h"

#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/Mu2eUtilities/inc/TriggerResultsNavigator.hh"

#include "TH1.h"
#include "TH2.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>

namespace mu2e {

  class CaloTriggerEffAna : public art::EDAnalyzer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag> caloClusterTag{Name("caloClusterTag"), Comment("Tag for CaloClusterCollection")};
      fhicl::Atom<art::InputTag> triggerResultsTag{Name("triggerResultsTag"), Comment("Tag for the global art::TriggerResults object")};
      fhicl::Atom<double> minEnergy{Name("minEnergy"), Comment("Minimum cluster energy for selection (MeV)")};
      fhicl::Atom<int> diagLevel{Name("diagLevel"), Comment("Diagnostic printouts"), 0};
    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit CaloTriggerEffAna(const Parameters& conf);

    void analyze(const art::Event& event) override;
    void beginJob() override;
    void endJob() override;

    void bookHistograms(const int index, std::string title);
    void fillHistograms(const int index, const CaloCluster& cluster);

  private:
    art::InputTag cluTag_;
    art::InputTag trigResultsTag_;
    double minEnergy_;
    int diagLevel_;

    const art::Event* event_;

    // Constant scaling factor for pileup rate calculation (1BB)
    const double kRateScaleFactor_ = 0.322 / 1.695e-6;

    // Global counters
    unsigned long long totalEvents_{0};
    unsigned long long totalOfflineEvents_{0};
    unsigned long long totalClustersPassingSelection_{0};

    // Path-specific map to store accepted counts
    std::map<std::string, unsigned long long> pathCounters_;
    std::map<std::string, unsigned long long> pathCountersOffline_;

    struct Hist_t {
      TH1F* cluster_energy{nullptr};
      TH1F* cluster_time{nullptr};
      TH2F* energy_vs_path{nullptr};
    };

    Hist_t hists_[100];
  };

  CaloTriggerEffAna::CaloTriggerEffAna(const Parameters& conf) :
    art::EDAnalyzer(conf),
    cluTag_(conf().caloClusterTag()),
    trigResultsTag_(conf().triggerResultsTag()),
    minEnergy_(conf().minEnergy()),
    diagLevel_(conf().diagLevel()) {}

  void CaloTriggerEffAna::beginJob() {
    bookHistograms(0, "all clusters");
    bookHistograms(1, "selected clusters");
  }

  void CaloTriggerEffAna::bookHistograms(const int index, std::string title) {
    Hist_t& Hist = hists_[index];
    art::ServiceHandle<art::TFileService> tfs;
    auto dir = tfs->mkdir(std::format("hist_{}", index), title);
    Hist.cluster_energy = dir.make<TH1F>("cluster_energy", "Deposited Energy;Energy [MeV];Clusters", 150, 0., 150.);
    Hist.cluster_time   = dir.make<TH1F>("cluster_time", "Time;Time [ns];Clusters", 200, 0., 2000.);
    Hist.energy_vs_path = dir.make<TH2F>("energy_vs_path", "Energy vs path;Path;Energy [MeV]", 50, 0, 50., 150, 0., 150.);
  }

  void CaloTriggerEffAna::fillHistograms(const int index, const CaloCluster& cluster) {
    Hist_t& Hist = hists_[index];
    if(!Hist.cluster_energy) {
      throw cet::exception("TRIGGER") << " Histogram set " << index << " not initialized!";
    }
    Hist.cluster_energy->Fill(cluster.energyDep());
    Hist.cluster_time->Fill(cluster.time());

    auto trigResultsHandle = event_->getValidHandle<art::TriggerResults>(trigResultsTag_);
    TriggerResultsNavigator nav(trigResultsHandle.product());
    const size_t npaths = nav.getTrigPaths().size();

    if(Hist.energy_vs_path->GetEntries() == 0) { // initialize the axes
      for(size_t i = 0; i < npaths; ++i) {
        const std::string path = nav.getTrigPathNameByIndex(i);
        Hist.energy_vs_path->Fill(path.c_str(), -1., 0.);
      }
    }

    for(size_t i = 0; i < npaths; ++i) {
      const std::string path = nav.getTrigPathNameByIndex(i);
      if (nav.accepted(path)) {
        Hist.energy_vs_path->Fill(path.c_str(), cluster.energyDep(), 1.);
      }
    }
  }

  void CaloTriggerEffAna::analyze(const art::Event& event) {
    totalEvents_++;
    event_ = &event;

    // Fetch Offline calorimeter clusters
    auto cluHandle = event.getValidHandle<CaloClusterCollection>(cluTag_);

    // Fetch the TriggerResults object
    auto trigResultsHandle = event.getValidHandle<art::TriggerResults>(trigResultsTag_);

    // Initialize the trigger navigator
    TriggerResultsNavigator nav(trigResultsHandle.product());
    if(diagLevel_ > 0) nav.print();

    // Track if the event passes the Offline selection
    bool pass = false;

    // Loop over clusters for counting and histogramming
    for (const auto& cluster : *cluHandle) {
      double cluEnergy = cluster.energyDep();

      // Log properties for all clusters before filtering
      fillHistograms(0, cluster);

      if (cluEnergy >= minEnergy_) {
        pass = true;
        totalClustersPassingSelection_++;
        fillHistograms(1, cluster);
      }
    }

    if(pass) ++totalOfflineEvents_;

    // Query trigger paths and count events passing the trigger
    const size_t npaths = nav.getTrigPaths().size();
    for (size_t i = 0; i < npaths; ++i) {
      const std::string path = nav.getTrigPathNameByIndex(i);
      if (nav.accepted(path)) {
        pathCounters_[path]++;
        if(diagLevel_ > 0) printf("  Accepting path %s\n", path.c_str());
        if(pass) pathCountersOffline_[path]++;
      }
    }
  }

  void CaloTriggerEffAna::endJob() {
    std::cout << "\n";
    std::cout << "=======================================================================================================\n";
    std::cout << "        Mu2e CaloTriggerEffAna Pileup Evaluation Report                  \n";
    std::cout << "=======================================================================================================\n";
    std::cout << " Total Events Processed: " << totalEvents_ << "\n";
    std::cout << " Total Clusters Passing Selection (E > " << minEnergy_ << " MeV): " << totalClustersPassingSelection_ << "\n";
    std::cout << " Rate Scaling Factor Applied: " << kRateScaleFactor_ << "\n";
    std::cout << "-------------------------------------------------------------------------------------------------------\n";
    std::cout << " " << std::left
              << std::setw(30) << "Path Name"
              << std::setw(15) << "N(Offline)"
              << std::setw(15) << "Eff(Offline)"
              << std::setw(15) << "All"
              << std::setw(15) << "Efficiency"
              << std::setw(12) << "Rate" << "\n";
    std::cout << "-------------------------------------------------------------------------------------------------------\n";

    for (const auto& [pathName, pathAccepts] : pathCounters_) {
      // Efficiency defined as N(events accepted) / N(events total)
      double efficiency = totalEvents_ > 0 ? static_cast<double>(pathAccepts) / totalEvents_ : 0.;
      auto nPassOffline = pathCountersOffline_.count(pathName) ? pathCountersOffline_[pathName] : 0;
      double eff_offline = totalOfflineEvents_ > 0 ? static_cast<double>(nPassOffline) / totalOfflineEvents_ : 0.;

      // Rate scaled by pileup factor
      double rate = efficiency * kRateScaleFactor_;

      std::cout << " " << std::left
                << std::setw(30) << pathName
                << std::setw(15) << nPassOffline
                << std::setw(15) << std::fixed << std::setprecision(5) << eff_offline
                << std::setw(15) << pathAccepts
                << std::setw(15) << std::fixed << std::setprecision(5) << efficiency
                << std::setw(12) << std::scientific << std::setprecision(4) << rate << "\n";
    }
    std::cout << "=======================================================================================================\n";
    std::cout << "\n";
  }
}

DEFINE_ART_MODULE(mu2e::CaloTriggerEffAna)
