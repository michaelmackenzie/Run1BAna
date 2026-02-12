// Dataset information
#ifndef __MU2E_EVT_ANA_DATASETS_C__
#define __MU2E_EVT_ANA_DATASETS_C__

namespace DATA {
  struct Dataset_t {
    Bool_t process_;
    TString name_;
    TString full_name_;
    Int_t    n_events_;
    Int_t    n_gen_events_;
  
    Dataset_t(Bool_t process, TString name, TString full_name, Int_t n_events, Int_t n_gen_events)
      : process_(process), name_(name), full_name_(full_name), n_events_(n_events), n_gen_events_(n_gen_events) {}
  };
  
  vector<Dataset_t> datasets() {
    vector<Dataset_t> datasets;
    datasets.emplace_back(true,  "cele0b1s5r0100", "nts.mu2e.CeEndpointMix1BBTriggered.MDC2025-000.root", 10000, 100000);

    return datasets;
  }
}
#endif // __MU2E_EVT_ANA_DATASETS_C__