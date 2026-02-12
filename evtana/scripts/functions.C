#ifndef __MU2E_EVT_ANA_FUNCTIONS_C__
#define __MU2E_EVT_ANA_FUNCTIONS_C__
#include "Mu2eEvtAna/scripts/datasets.C"

Mu2eEvtAna::Mu2eEvtAna* gMu2eAna = nullptr;
Mu2eEvtAna::RMCAna* gRMCAna = nullptr;

// Using the default Mu2eEvtAna analyzer
int mu2e_ana(TString dataset, int Mode = 0, Long64_t max_entries = 1e6, Long64_t first_entry = 0) {

  // Find the dataset
  TString file_list = "";
  for(auto config : DATA::datasets()) {
    if(config.name_ == dataset) {
      file_list = Form("Mu2eEvtAna/file_lists/%s.files", config.full_name_.Data());
    }
  }
  if(file_list == "") {
    cout << "Dataset " << dataset << " not found!" << endl;
    return -1;
  }

  // Setup a Mu2eEvtAna processing
  if(gMu2eAna) delete gMu2eAna;
  gMu2eAna = new Mu2eEvtAna::Mu2eEvtAna(0);
  gMu2eAna->AddFile(file_list, max_entries, first_entry);
  gMu2eAna->SetName(Form("mu2e_ana_%s_m%i", dataset.Data(), Mode));
  gMu2eAna->cache_size_ = 200000000U; //200 MB cache
  // gMu2eAna->cache_size_ = 2000000000U; //2 GB cache
  gMu2eAna->load_baskets_ = false;
  gMu2eAna->report_rate_ = 5000;

  // Process the tree
  const int status = gMu2eAna->Process(max_entries);
  cout << "Status code = " << status << endl;

  return status;
}

// Using the RMCAna analyzer
int rmc_ana(TString dataset, int Mode = 0, Long64_t max_entries = -1, Long64_t first_entry = 0) {

  // Find the dataset
  TString file_list = "";
  for(auto config : DATA::datasets()) {
    if(config.name_ == dataset) {
      file_list = Form("Mu2eEvtAna/file_lists/%s.files", config.full_name_.Data());
    }
  }
  if(file_list == "") {
    cout << "Dataset " << dataset << " not found!" << endl;
    return -1;
  }

// Setup a RMCAna processing
if(gRMCAna) delete gRMCAna;
gRMCAna = new Mu2eEvtAna::RMCAna(0);
gRMCAna->AddFile(file_list, max_entries, first_entry);
gRMCAna->SetName("rmc_ana");
gRMCAna->cache_size_ = 2000000000U; //2 GB cache

// Process the tree
const int status = gRMCAna->Process(max_entries);
  cout << "Status code = " << status << endl;

  return status;
}

#endif // __MU2E_EVT_ANA_FUNCTIONS_C__