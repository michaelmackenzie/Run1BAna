#include "Run1BAna/evtana/inc/ConvAna.hh"

using namespace mu2e;
namespace Run1BEvtAna {

  //------------------------------------------------------------------------------------
  // Constructor
  ConvAna::ConvAna(int verbose) : Run1BEvtAna(verbose) {
    //initialize the arrays to null
    for(int ihist = 0; ihist < kMaxHists; ++ihist) {
      sys_hists_[ihist] = nullptr;
    }
  }


  //------------------------------------------------------------------------------------
  // Define the histogram selections
  void ConvAna::InitHistSelections() {

    //-----------------------------------------------------------------------------
    // book histogram selections
    //-----------------------------------------------------------------------------
    struct hist_info_t {
      TString _dsc; // description of the selection
      bool    _trk; // track histograms
      bool    _hlx; // helix histograms
      bool    _smp; // sim particle histograms
      bool    _gnp; // gen particle histograms
      bool    _crv; // CRV histograms
      bool    _sys; // systematic histograms
      bool    _crs; // control regions included
      bool    _trs; // output trees
      hist_info_t(TString dsc = "", bool trk = false, bool hlx = false, bool smp = false, bool gnp = false,
                  bool crv = false, bool sys = false, bool crs = false, bool trs = false)
        : _dsc(dsc), _trk(trk), _hlx(hlx), _smp(smp), _gnp(gnp), _crv(crv), _sys(sys),
          _crs(crs), _trs(trs) {}
    };

    hist_info_t* hist_sets[kMaxHists];
    for (int i=0; i<kMaxHists; i++) {
      hist_sets[i] = nullptr;
    }

    //                                 description                         trk    hlx    simp   genp   crv    sys    crs    trs
    hist_sets[  0] = new hist_info_t("All events, Offline track"        ,  true,  true,  true,  true,  true, false, false, false);
    hist_sets[  1] = new hist_info_t("All events, APR track"            ,  true,  true,  true,  true, false, false, false, false);
    hist_sets[  2] = new hist_info_t("All events, CPR track"            ,  true,  true,  true,  true, false, false, false, false);
    hist_sets[  3] = new hist_info_t("No weights, Offline track"        ,  true,  true,  true,  true, false, false, false, false);
    hist_sets[  4] = new hist_info_t("e- and p > 95"                    ,  true, false,  true, false, false, false, false, false);
    hist_sets[  5] = new hist_info_t("e+ and p > 75"                    ,  true, false,  true, false, false, false, false, false);
    hist_sets[  6] = new hist_info_t("Gen(E) > 95"                      ,  true, false,  true, false, false, false, false, false);
    hist_sets[  7] = new hist_info_t("e+-, TrkID, no MC cut"            ,  true, false,  true, false, false, false, false, false);
    hist_sets[  8] = new hist_info_t("Events with a track"              ,  true, false,  true, false, false, false, false, false);
    hist_sets[ 10] = new hist_info_t("Offline, event selection"         ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 11] = new hist_info_t("e+/-: trigger"                    ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 12] = new hist_info_t("e+/-: failed trigger"             ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 13] = new hist_info_t("e+/-: negative"                   ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 14] = new hist_info_t("e+/-: positive"                   ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 15] = new hist_info_t("e+/-: no weights"                 ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 16] = new hist_info_t("e-: p > 95"                       ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 20] = new hist_info_t("e-: full window"                  ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 21] = new hist_info_t("e-: narrow window"                ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 22] = new hist_info_t("e-: full window, no weights"      ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 23] = new hist_info_t("e-: high error"                   ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 25] = new hist_info_t("e-: mom up"                       ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 26] = new hist_info_t("e-: mom down"                     ,  true, false, false, false, false, false,  true, false);
    hist_sets[ 30] = new hist_info_t("e-: no CRV veto"                  ,  true, false, false, false,  true, false,  true, false);
    hist_sets[ 31] = new hist_info_t("e+-: no CRV veto or p cut"        ,  true, false, false, false,  true, false,  true, false);
    hist_sets[ 34] = new hist_info_t("e-: low dP(ST)"                   ,  true, false, false, false, false,  true,  true, false);
    hist_sets[ 35] = new hist_info_t("e-: high dP(ST)"                  ,  true, false, false, false, false,  true,  true, false);
    hist_sets[ 40] = new hist_info_t("e+: full window"                  ,  true,  true,  true,  true,  true,  true,  true,  true);
    hist_sets[ 41] = new hist_info_t("e+: narrow window"                ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 42] = new hist_info_t("e+: full window, no weights"      ,  true,  true,  true,  true, false, false,  true, false);
    hist_sets[ 50] = new hist_info_t("e+: no CRV veto"                  ,  true, false, false, false,  true, false,  true, false);

    // CRV studies histograms
    hist_sets[ 80] = new hist_info_t("CRV: 1"                           ,  true, false, false, false,  true, false, false, false);
    hist_sets[ 81] = new hist_info_t("CRV: Ue/Umu tracks"               ,  true, false, false, false,  true, false, false, false);
    hist_sets[ 82] = new hist_info_t("CRV: calo cluster with track"     ,  true, false, false, false,  true, false, false, false);
    hist_sets[ 83] = new hist_info_t("CRV: dt < -60 ns"                 ,  true, false, false, false,  true, false, false, false);
    hist_sets[ 84] = new hist_info_t("CRV: MC electrons"                ,  true, false, false, false,  true, false, false,  true);
    hist_sets[ 85] = new hist_info_t("CRV: MC muons"                    ,  true, false, false, false,  true, false, false,  true);
    hist_sets[ 86] = new hist_info_t("CRV: Correct cluster"             ,  true, false, false, false,  true, false, false, false);
    hist_sets[ 87] = new hist_info_t("CRV: Upstream veto"               ,  true, false, false, false,  true, false, false, false);

    for (int i=0; i<kMaxHists; i++) {
      const int index = i % 1000;
      if(!hist_sets[index]) continue;
      const bool is_cr = i >= 1000 && !hist_sets[index]->_crs;
      if(is_cr) continue; // control region
      if(i >= 4000) break; //Control regions above 4000 not yet implemented
      evt_hists_[index] = new EventHist_t;
      if(hist_sets[index]->_trk) trk_hists_[index] = new TrackHist_t;
      if(hist_sets[index]->_crv && ! is_cr) crv_hists_[index] = new CRVHist_t;
      if(hist_sets[index]->_sys) sys_hists_[index] = new SysHist_t;
      // FIXME: Add missing histogram types
    }
  }

  //------------------------------------------------------------------------------------
  // Initialize systematic histograms
  void ConvAna::BookSystematicHist(SysHist_t* Hist, const char* Folder) {
    if(!Hist) {
      throw std::runtime_error("Attempting to book histograms in a null SysHist_t\n");
    }
    for(int isys = 0; isys < kMaxSystematics; ++isys) {
      // check if the systematic is defined
      TString name = systematics_.GetName(isys);
      if(name == "") continue;
      Hist->fObs[isys] = new TH1F(Form("obs_%i", isys),Form("%s: Systematic %s",Folder, name.Data()), 150, 80., 110.); //FIXME: This should inherit from the nominal observable binning
      // For debug investigations
      if(fill_verbose_sys_) {
        Hist->fDeltaObs   [isys] = new TH1F(Form("delta_obs_%i"   , isys),Form("%s: Systematic %s: #DeltaObs"   ,Folder, name.Data()), 100, -2., 2.);
        Hist->fWeight     [isys] = new TH1F(Form("weight_%i"      , isys),Form("%s: Systematic %s: Weight"      ,Folder, name.Data()), 100,  0., 2.);
        Hist->fDeltaWeight[isys] = new TH1F(Form("delta_weight_%i", isys),Form("%s: Systematic %s: #DeltaWt/Wt" ,Folder, name.Data()), 100, -1., 1.);
      }
    }
  }

  //------------------------------------------------------------------------------------
  // Initialize the histogram sets
  void ConvAna::BookHistograms(TDirectory* dir) {
    Run1BEvtAna::BookHistograms(dir);
    for(int ihist = 0; ihist < kMaxHists; ++ihist) {
      if(sys_hists_[ihist]) {
        const char* folder = Form("sys_%i", ihist);
        auto subdir = dir->mkdir(folder);
        subdir->cd();
        BookSystematicHist(sys_hists_[ihist], folder);
        dir->cd();
        sys_dirs_[ihist] = subdir;
      }
    }
  }

  //------------------------------------------------------------------------------------
  // Initialize the input ntuple information
  int ConvAna::InitializeInput() {
    Run1BEvtAna::InitializeInput();
    return 0;
  }

  //------------------------------------------------------------------------------------
  // Initialize the output ntuple information
  int ConvAna::InitializeOutput() {
    Run1BEvtAna::InitializeOutput();
    return 0;
  }


  //------------------------------------------------------------------------------------
  // Initialize event information
  void ConvAna::InitializeEvent() {
    Run1BEvtAna::InitializeEvent();
  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool ConvAna::ProcessEvent() {
    FillEventHist(evt_hists_[0]); //all events with well defined inputs
    if(evt_.nde_tracks_ != 1) return false; //exactly one positron or electron
    FillEventHist(evt_hists_[1]);
    return true;
  }
}
