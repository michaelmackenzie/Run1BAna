// Test the RMCAna selector

using namespace mu2e::Mu2eEvtAna;

int test_rmc(Long64_t max_entries = 10000) {

  // Retrieve an EventNtuple tree
  auto f_in = TFile::Open("/exp/mu2e/app/users/mmackenz/Yale/main/nts.owner.description.version.sequencer.root", "READ");
  if(!f_in) return -1;
  auto t_in = (TTree*) f_in->Get("EventNtuple/ntuple");
  if(!t_in) {
    cout << "Input ntuple not found!\n";
    f_in->ls();
    return 1;
  }

  // Setup a Mu2eEvtAna processing
  RMCAna ana(0);
  ana.SetInput(t_in);
  ana.SetName("test");

  // Process the tree
  const int status = ana.Process(max_entries);
  cout << "Status code = " << status << endl;

  return 0;
}
