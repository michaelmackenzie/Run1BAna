// Train a photon vs. pileup model

// Make a plot of response vs. variable
int response_plot(TTree* tree, TString var, TString mva, TString name, bool signal = true) {
  TCanvas c;
  tree->Draw(mva+":"+var, (signal) ? "classID == 0" : "classID == 1", "colz");
  c.SaveAs(name + ".png");
  c.SetLogz();
  c.SaveAs(name + "_log.png");
  return 0;
}

// Main function
int train_mva(bool train = true) {

  // Get the input files
  TFile* f_sig = TFile::Open("/exp/mu2e/data/users/mmackenz/run1b/histograms/nts.mu2e.FlatGammaMix1BB-KL.Run1Bah_best_v1_4-001.root", "READ");
  TFile* f_bkg = TFile::Open("/exp/mu2e/data/users/mmackenz/run1b/histograms/nts.mu2e.NoPrimaryMix1BB-KL_skim_clusters.Run1Bah_best_v1_4-001.root", "READ");
  if(!f_bkg || !f_sig) return 1;

  // Retrieve the input trees
  TTree* t_sig = (TTree*) f_sig->Get("Run1BAna/tree_60/tree"); // Base photon selection
  TTree* t_bkg = (TTree*) f_bkg->Get("Run1BAna/tree_60/tree");
  if(!t_sig || !t_bkg) {
    cout << "Input trees not found!\n";
    return 2;
  }

  /////////////////////////////////////////////////////////////////////////
  // ---------------- Now we use ROOT TMVA to train -------------------- //
  /////////////////////////////////////////////////////////////////////////

  TMVA::Tools::Instance();

  TString tname = "train_mva_output"; // training name
  gSystem->Exec(Form("[ ! -d %s ] && mkdir %s", tname.Data(), tname.Data()));
  gSystem->cd(tname.Data());

  TString outfileName = tname + ".root";
  TString outFolder = "tmva_" + tname;

  if(train) {
    TFile* outfile = new TFile(outfileName.Data(), "RECREATE");

    // instantiate TMVA::DataLoader object
    TMVA::DataLoader* dataLoader = new TMVA::DataLoader(tname);

    // instantiate TMVA::Factory object
    TMVA::Factory* factory = new TMVA::Factory("TMVAClassification", outfile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

    // define variables to be used for training
    dataLoader->AddVariable ("cluster_frac_1"        , "E_{1}/E"      , ""       , 'F');
    dataLoader->AddVariable ("cluster_frac_2"        , "E_{1+2}/E"    , ""       , 'F');
    dataLoader->AddVariable ("cluster_t_var"         , "#sigma^{2}(t)", "ns^{2}" , 'F');
    dataLoader->AddVariable ("cluster_second_moment" , "Second moment", "MeV"    , 'F');
    dataLoader->AddVariable ("cluster_radius"        , "Radius"       , "mm"     , 'F');
    dataLoader->AddVariable ("cluster_ncr"           , "N(crystals)"  , ""       , 'F');

    dataLoader->AddSpectator("cluster_energy"        , "Energy"       , "MeV"    , 'F');
    dataLoader->AddSpectator("cluster_time"          , "Time"         , "ns"     , 'F');
    dataLoader->AddSpectator("event_weight"          , "Weight"       , ""       , 'F');
    dataLoader->AddSpectator("gen_energy"            , "Gen energy"   , "MeV"    , 'F');

    // register the trees
    dataLoader->AddSignalTree    (t_sig, 1.);
    dataLoader->AddBackgroundTree(t_bkg, 1.);

    // set training options
    // setting all to zero allows half of tracks to be used for training, and half for testing
    const double training_fraction = 0.7; // use 70% of the data for training
    const double max_train_events  = 100000.; // no need to use more than 100k training events
    const int nsignal     = std::min(max_train_events, t_sig->GetEntries() * training_fraction);
    const int nbackground = std::min(max_train_events, t_bkg->GetEntries() * training_fraction);
    TString trainingOptions = Form("nTrain_Signal=%i:nTrain_Background=%i", nsignal, nbackground);
    trainingOptions += ":nTest_Signal=0:nTest_Background=0";
    trainingOptions += ":SplitMode=Random:!V";

    // final preparation of training and test trees
    TCut preselectionCut = ""; // empty since preselection is performed in created sig and bkg trees
    dataLoader->PrepareTrainingAndTestTree(preselectionCut, trainingOptions);

    // book MVA methods
    factory->BookMethod(dataLoader, TMVA::Types::kMLP, "MLP", "!H:!V:VarTransform=N:NeuronType=tanh:HiddenLayers=N+1,N:TestRate=5:!UseRegulator");
    factory->BookMethod(dataLoader, TMVA::Types::kBDT, "BDT",
                        "!H:!V:NTrees=850:nEventsMin=150:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning");

    // now we perform training, testing, and evaluation of all methods booked
    factory->SetVerbose();
    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();

    // close output file
    outfile->Close();

    // cleanup
    delete factory;
    delete dataLoader;
  }

  // make additional figures
  if(!gROOT->IsBatch()) TMVA::TMVAGui(outfileName); // if interactive, use the GUI
  else {
    TMVA::variables   (tname.Data(), outfileName.Data(), "InputVariables_Id"); //1D variable plots
    TMVA::mvas        (tname.Data(), outfileName.Data(), TMVA::kCompareType); //Training + Testing MVA scores
    TMVA::correlations(tname.Data(), outfileName.Data()); //linear correlations
    TMVA::efficiencies(tname.Data(), outfileName.Data()); //ROC

    TFile* f = TFile::Open(outfileName, "READ");
    if(!f) return 1;
    TTree* t_test = (TTree*) f->Get(tname + "/TestTree");
    if(!t_test) {f->Close(); return 2;}
    TString figdir = tname + "/plots/";
    response_plot(t_test, "cluster_energy", "MLP", figdir + "energy_vs_MLP_test");
    f->Close();
  }

  return 0;

}
