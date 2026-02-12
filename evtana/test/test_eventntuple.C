// Test Mu2e/EventNtuple format

#include "EventNtuple/inc/TrkInfo.hh"

#include "TFile.h"
#include "TTree.h"

int test_eventntuple(const Long64_t nprocess = 1) {

  // Retrieve an EventNtuple tree
  auto f_in = TFile::Open("/pnfs/mu2e/tape/phy-nts/nts/mu2e/CeEndpointOnSpillTriggered/MDC2020aq_best_v1_3_v06_03_00/root/ed/df/nts.mu2e.CeEndpointOnSpillTriggered.MDC2020aq_best_v1_3_v06_03_00.001210_00000000.root", "READ");
  if(!f_in) return -1;
  auto t_in = (TTree*) f_in->Get("EventNtuple/ntuple");
  if(!t_in) {
    cout << "Input ntuple not found!\n";
    f_in->ls();
    return 1;
  }
  const Long64_t nentries = t_in->GetEntries();
  const Long64_t max_entry = min(nprocess, nentries); //only process the requested events

  // Test reading a list of tracks in each event
  vector<mu2e::TrkInfo>* trks = nullptr;
  vector<vector<mu2e::TrkSegInfo>>* trksegs = nullptr;
  t_in->SetBranchAddress("trk", &trks);
  t_in->SetBranchAddress("trksegs", &trksegs);

  // Test reading the track quality
  vector<mu2e::MVAResultInfo>* trk_qual = nullptr;
  t_in->SetBranchAddress("trkqual", &trk_qual);

  // Test reading the hit counts
  mu2e::HitCount* hit_cnt = nullptr;
  t_in->SetBranchAddress("hitcount.", &hit_cnt);

  // Process the requested events
  for(int entry = 0; entry < max_entry; ++entry) {
    cout << "Beginning to process entry " << entry << endl;
    if(!t_in) {
      cout << "  Tree is no longer defined!" << endl;
      break;
    }
    t_in->GetEntry(entry);
    if(trks) {
      cout << "  Tracks not null, size = " << trks->size() << endl;
      for(auto trk : *trks) {
        cout << "  --> Track: status = " << trk.status << ", pdg = " << trk.pdg << endl;
      }
    }
    if(trk_qual) {
      cout << "  Track quality not null, size = " << trk_qual->size() << endl;
      for(auto qual : *trk_qual) {
        cout << "  --> Score: valid = " << qual.valid << ", score = " << qual.result << endl;
      }
    }
    if(trksegs) {
      cout << "  Track segments not null, size = " << trksegs->size() << endl;
      for(auto trkseg : *trksegs) {
        cout << "  --> Track segments:" << endl;
        for(auto seg : trkseg) {
          cout << "    --> Track: p = " << seg.mom.r() << ", z = " << seg.pos.z() << endl;
        }
      }
    }
    if(hit_cnt) {
      cout << "  Hit information:" << endl;
      cout << "    Digis: " << hit_cnt->nsd << endl;
    }
  }

  return 0;
}
