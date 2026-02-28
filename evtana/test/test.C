// Test the Run1BEvtAna tools
Run1BEvtAna::Run1BEvtAna* gAna;

int test(Long64_t max_entries = 1e4, Long64_t first_entry = 0) {

  // TString file_list = "/pnfs/mu2e/tape/phy-nts/nts/mu2e/CeEndpointOnSpillTriggered/MDC2020aq_best_v1_3_v06_03_00/root/ed/df/nts.mu2e.CeEndpointOnSpillTriggered.MDC2020aq_best_v1_3_v06_03_00.001210_00000000.root";
  TString file_list = "Run1BAna/file_lists/nts.mu2e.CeEndpointMixLowTriggerable-KL.Run1B-003.root.files";

  // Setup a Mu2eEvtAna processing
  gAna = new Run1BEvtAna::Run1BEvtAna((max_entries < 10) ? (max_entries == 1) ? 10 : 2 : 0);
  gAna->AddFile(file_list, max_entries, first_entry);
  gAna->SetName("test");
  gAna->load_baskets_ = false;
  gAna->cache_size_ = 2000000000U; //2 GB cache

  // Process the tree
  const int status = gAna->Process(max_entries);
  cout << "Status code = " << status << endl;

  return 0;
}
