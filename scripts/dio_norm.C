// Get the DIO normalization information

//---------------------------------------------------------------------------
// DIO theoretical spectrum
TH1* get_dio_spectrum() {
  TTree tree("t1","t1");
  TString table = "Run1BAna/data/heeck_finer_binning_2016_szafron.tbl";
  int     nb    = 11000; //finer binning in this table
  double  bin   = 0.01;
  tree.ReadFile(table,"e/D:w/D");
  int n = tree.GetEntries();

  double emin = 0.;
  double emax = 110.;
  TH1* h_dio = new TH1D("h_dio","DIO spectrum",nb,emin,emax);

  double e, w;

  tree.SetBranchAddress("e",&e);
  tree.SetBranchAddress("w",&w);

  int prev_bin = 0;
  for (int i = 0; i < n; ++i) {
    tree.GetEntry(i);
    int ibin = h_dio->FindBin(e-bin/2.);
    if(prev_bin > 0 && prev_bin != ibin - 1) {
      cout << "Bin " << ibin << " entry " << i << " E = " << e
           << " but prev_bin = " << prev_bin << endl;
    }
    h_dio->SetBinContent(ibin,w);
    prev_bin = ibin;
  }
  h_dio->SetLineColor(kBlue);
  h_dio->SetLineWidth(2);
  h_dio->SetFillColor(0);
  h_dio->SetFillStyle(0);
  h_dio->Scale(1./(bin*h_dio->Integral()));
  return h_dio;
}

//---------------------------------------------------------------------------
// evaluate the integral in a range
double dio_norm_range(TH1* h_dio, double emin, double emax) {
  return h_dio->Integral(h_dio->FindBin(emin), h_dio->FindBin(emax))*h_dio->GetBinWidth(1);
}

//---------------------------------------------------------------------------
// main function
void dio_norm() {
  // Get the PDFs
  TH1* h_dio = get_dio_spectrum();

  const double dio_0_inf  = dio_norm_range(h_dio,  0., 105.);
  const double dio_0_60   = dio_norm_range(h_dio,  0.,  60.);
  const double dio_60_80  = dio_norm_range(h_dio, 60.,  80.);
  const double dio_80_90  = dio_norm_range(h_dio, 80.,  90.);
  const double dio_60_inf = dio_norm_range(h_dio, 60., 105.);
  const double dio_80_inf = dio_norm_range(h_dio, 80., 105.);
  const double dio_90_inf = dio_norm_range(h_dio, 90., 105.);
  const double dio_95_inf = dio_norm_range(h_dio, 95., 105.);

  // Print out the integrals
  cout << "-----------------------------------------------------------------\n"
       << "DIO  0 - inf  = " << dio_0_inf  << endl
       << "DIO  0 - 60   = " << dio_0_60   << endl
       << "DIO 60 - 80   = " << dio_60_80  << endl
       << "DIO 80 - 90   = " << dio_80_90  << endl
       << "DIO 60 - inf  = " << dio_60_inf << endl
       << "DIO 70 - inf  = " << dio_80_inf << endl
       << "DIO 90 - inf  = " << dio_90_inf << endl
       << "DIO 95 - inf  = " << dio_95_inf << endl
       << "-----------------------------------------------------------------\n";
}
