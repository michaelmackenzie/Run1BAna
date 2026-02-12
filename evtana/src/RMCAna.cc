#include "Run1BAna/evtana/inc/RMCAna.hh"

using namespace mu2e;
namespace Run1BEvtAna {

  //------------------------------------------------------------------------------------
  // Constructor
  RMCAna::RMCAna(int verbose) : Run1BEvtAna(verbose) {
  }


  //------------------------------------------------------------------------------------
  // Define the histogram selections
  void RMCAna::InitHistSelections() {
    //Default histogram selections
    evt_hists_[0] = new EventHist_t;
    trk_hists_[0] = new TrackHist_t;
    crv_hists_[0] = new CRVHist_t;

    // Basic RMC selection
    evt_hists_[1] = new EventHist_t;
    trk_hists_[1] = new TrackHist_t;
  }

  //------------------------------------------------------------------------------------
  // Initialize the input ntuple information
  int RMCAna::InitializeInput() {
    Run1BEvtAna::InitializeInput();
    return 0;
  }

  //------------------------------------------------------------------------------------
  // Initialize the output ntuple information
  int RMCAna::InitializeOutput() {
    Run1BEvtAna::InitializeOutput();
    return 0;
  }


  //------------------------------------------------------------------------------------
  // Initialize event information
  void RMCAna::InitializeEvent() {
    Run1BEvtAna::InitializeEvent();
  }

  //------------------------------------------------------------------------------------
  // Main event-by-event processing
  bool RMCAna::ProcessEvent() {
    FillEventHist(evt_hists_[0]); //all events with well defined inputs
    if(evt_.nde_tracks_ != 1) return false; //exactly one positron or electron
    FillEventHist(evt_hists_[1]);
    return true;
  }
}
