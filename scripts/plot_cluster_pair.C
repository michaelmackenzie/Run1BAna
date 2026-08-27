// Plot cluster energy from first + second disk
#include "EventNtuple/inc/CaloClusterInfo.hh"
using namespace mu2e;

int plot_cluster_pair(TString file_name) {
  TFile* f = TFile::Open(file_name, "READ");
  if(!f) return 1;

  TTree* t = (TTree*) f->Get("EventNtuple/ntuple");
  if(!t) {
    cout << "No ntuple found!\n";
    f->ls();
    return 1;
  }

// *Br  331 :caloclusters.time_ : Float_t time_[caloclusters_]                  *
// *Entries :    47001 : Total  Size=     360837 bytes  File Size  =     216030 *
// *Baskets :       10 : Basket Size=     393216 bytes  Compression=   1.67     *
// *............................................................................*
// *Br  332 :caloclusters.timeErr_ : Float_t timeErr_[caloclusters_]            *
// *Entries :    47001 : Total  Size=     360879 bytes  File Size  =     215501 *
// *Baskets :       10 : Basket Size=     393216 bytes  Compression=   1.67     *
// *............................................................................*
// *Br  333 :caloclusters.energyDep_ : Float_t energyDep_[caloclusters_]        *
// *Entries :    47001 : Total  Size=     360907 bytes  File Size  =     212939 *
// *Baskets :       10 : Basket Size=     393216 bytes  Compression=   1.69     *
// *............................................................................*
// *Br  334 :caloclusters.energyDepErr_ : Float_t energyDepErr_[caloclusters_]  *
// *Entries :    47001 : Total  Size=     360949 bytes  File Size  =     209906 *
// *Baskets :       10 : Basket Size=     393216 bytes  Compression=   1.72     *
// *............................................................................*
// *Br  335 :caloclusters.cog_.fCoordinates.fX : Float_t fX[caloclusters_]      *
// *Entries :    47001 : Total  Size=     361011 bytes  File Size  =     222223 *
// *Baskets :       10 : Basket Size=     393216 bytes  Compression=   1.62     *
// *............................................................................*
// *Br  336 :caloclusters.cog_.fCoordinates.fY : Float_t fY[caloclusters_]      *
// *Entries :    47001 : Total  Size=     361011 bytes  File Size  =     215523 *
// *Baskets :       10 : Basket Size=     393216 bytes  Compression=   1.67     *
// *............................................................................*
// *Br  337 :caloclusters.cog_.fCoordinates.fZ : Float_t fZ[caloclusters_]      *
  std::vector<CaloClusterInfo>* clusters = new std::vector<CaloClusterInfo>;
  t->SetBranchAddress("caloclusters"     , &clusters);
  for(Long64_t entry = 0; entry < t->GetEntriesFast(); ++entry) {
    t->GetEntry(entry);
    const auto nclusters = clusters->size();
    if(nclusters > 1) {
      cout << "Entry " << entry << ": " << nclusters << " clusters\n";
      for(size_t icl = 0; icl < nclusters; ++icl) {
        const auto& cluster = clusters->at(icl);
        printf("  %2zu: %i %.1f %.1f %.1f %.1f\n", icl, cluster.diskID_, cluster.energyDep_, cluster.time_, cluster.cog_.X(), cluster.cog_.Y());
      }
    }
    if(entry >= 1000) break;
  }

  return 0;
}
