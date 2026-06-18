// Plot muon distributions in the virtual detectors
void plotMuonsVD(const bool ds_on = true) {
  TString ds_off_path = "/scratch/mu2e/users/jmott/PlanB/CaloStopsTS5/ExtraVDs/FieldOff/outputTrees_Run1B_CaloStops.root";
  TString ds_on_path  = "/scratch/mu2e/users/jmott/PlanB/CaloStopsTS5/ExtraVDs/FieldOn/outputTrees_Run1A_TargetStops.root";

  TString file_name = (ds_on) ? ds_on_path : ds_off_path;

  // Fixed number of virtual detector planes
  const int nVDs = 14; // There are 14 planes here - 11 DS entrance planes, 2 Col5 planes, but only 1 ST planes
  const int nMaxVDs = 15; // This how big it could be from the filling of arrays

  // Setup to Read data from file
  TFile* file = TFile::Open(file_name, "READ");
  if(!file) return;
  TTree* tree = (TTree*) file->Get("stepPointMCDumperMuonVDs/nt");
  if(!tree) return;

  // Set branch structure
  unsigned nVDs_;
  unsigned vdID_[nMaxVDs];
  float x_[nMaxVDs];
  float y_[nMaxVDs];
  float z_[nMaxVDs];
  float px_[nMaxVDs];
  float py_[nMaxVDs];
  float pz_[nMaxVDs];
  float t_[nMaxVDs];
  tree->SetBranchAddress("nVDs",&nVDs_);
  tree->SetBranchAddress("vdID",&vdID_);
  tree->SetBranchAddress("x",x_);
  tree->SetBranchAddress("y",y_);
  tree->SetBranchAddress("z",z_);
  tree->SetBranchAddress("px",px_);
  tree->SetBranchAddress("py",py_);
  tree->SetBranchAddress("pz",pz_);
  tree->SetBranchAddress("t",t_);

  // Histograms
  TH2F* hXY[nMaxVDs];
  for(unsigned i = 0; i < nVDs; ++i) {
    hXY[i] = new TH2F(Form("xy_%i", i), "Muon (x, y)", 100, -4200., -3700., 100, -200., 200.);
  }

  // Loop over file and fill histograms
  int entries = tree->GetEntries();
  float target = 0.0;
  bool first = true;
  for(int iEvt = 0; iEvt < entries; iEvt++){
    if(iEvt % 10000 == 0) cout << "Processing event " << iEvt << "...\n";

    // Read next event
    tree->GetEntry(iEvt);

    // Check we had hits in all planes
    if(nVDs_ != nVDs) continue;

    if(first) {
      for(unsigned i = 0; i < nVDs; ++i) {
        hXY[i]->SetTitle(Form("Muon (x, y) - z = %.0f mm;x (mm); y (mm)", z_[i]));
      }
      first = false;
    }

    // Loop over all VDs
    for(unsigned i = 0; i < nVDs; ++i) {
      hXY[i]->Fill(x_[i], y_[i]);
    }

  }

  for(unsigned i = 0; i < nVDs; ++i) {
    new TCanvas();
    hXY[i]->Draw("colz");
  }
}
