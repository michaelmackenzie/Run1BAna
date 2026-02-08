#ifndef __RUN1BANA_ANALYSIS_DATASETINFO__
#define __RUN1BANA_ANALYSIS_DATASETINFO__
#include "Run1BAna/analysis/types.C"
#include "Run1BAna/analysis/physics.C"
#include <map>

std::map<TString, DatasetInfo_t> datasets_;
const char* hist_func_ = "run1b_ana"; // for histogram file naming
int         hist_mode_ = 2;

DatasetInfo_t get_dataset_info(TString name) {
  // Initialize the datasets if not already
  if(datasets_.size() == 0) {
    const double rate_dio((1.-muon_capture_fraction_)*nmuons_per_pot_), rate_dio_50(rate_dio*dio_frac_50_),
      rate_dio_80(rate_dio*dio_frac_80_), rate_dio_90(rate_dio*dio_frac_90_), rate_dio_95(rate_dio*dio_frac_95_);
    const double rate_rpc(pion_stop_rate_*pion_survive_frac_*rpc_br_), rate_rpc_int(rate_rpc*rpc_int_br_),
      rate_rmc(nmuons_per_pot_*muon_capture_fraction_*br_rmc_/rmc_frac_57_), rate_rmc_conv(rate_rmc*rmc_conv_),
      rate_rmc_int(rate_rmc*rmc_int_br_); // FIXME: rmc_int_br in convolution sometimes
    const double rate_ipa_dio((1.-muon_capture_carbon_)*ipa_nmuons_per_pot_*ipa_dio_frac_70_);
    const double rate_rmc_85(rate_rmc*rmc_frac_85_), rate_rmc_85_int(rate_rmc_85*rmc_int_br_);
    const double rate_pbar(pbar_stops_per_pot_);
    const double rate_pileup(1.);
    const double rate_sig(muon_capture_fraction_*nmuons_per_pot_);
    // Retrieve N(gen) using scripts/samCountGenEvents.sh and N(events) using scripts/samCountEvents.sh
    //                                  N(gen)      N(digi) emin emax    rate          stn dataset           art dataset

    datasets_["dio_50"]  = DatasetInfo_t(  1e6 ,   9172, 50., 105., rate_dio_50   , "dio50b1s5r0100", "nts.mu2e.DIOtail50MixLowTriggerable-KL.Run1Bab2_best_v1_2.root");
    datasets_["dio_80"]  = DatasetInfo_t(  1e7 , 223629, 80., 105., rate_dio_80   , "dio80b1s5r0100", "nts.mu2e.DIOtail80MixLowTriggerable-KL.Run1Bab2_best_v1_2.root");
    datasets_["dio_90"]  = DatasetInfo_t( 9.3e6, 225190, 90., 105., rate_dio_90   , "dio90b1s5r0100", "nts.mu2e.DIOtail90MixLowTriggerable-KL.Run1Bab2_best_v1_2.root");
    datasets_["pileup"]  = DatasetInfo_t(  1e6 , 204443,  0.,   1., rate_pileup   , "mnbs0b1s5r0100", "nts.mu2e.NoPrimaryMixLowTriggerable-KL.Run1Bab2_best_v1_2.root");
    // datasets_["mumem"]   = DatasetInfo_t( 9.3e6, 225190, 90., 105., rate_sig      , "dio90b1s5r0100", "nts.mu2e.DIOtail90MixLowTriggerable-KL.Run1Bab2_best_v1_2.root");
    // datasets_["mumem"]   = DatasetInfo_t( 1e7 , 223629, 80., 105., rate_sig      , "dio80b1s5r0100", "nts.mu2e.DIOtail80MixLowTriggerable-KL.Run1Bab2_best_v1_2.root");
    datasets_["mumem"]   = DatasetInfo_t( 1e6 ,  22717, 0.,   1., rate_sig       , "cele0b1s5r0100", "nts.mu2e.CeEndpointMixLowTriggerable-KL.Run1Bab2_best_v1_2.root");

  }
  if(datasets_.count(name) != 0) return datasets_[name];
  cout << __func__ << ": No dataset with name " << name.Data() << " found!\n";
  return DatasetInfo_t();
}

//---------------------------------------------------------------------------------------------------------------------------
void set_style(const TString name, TString& title, int& color) {
  title = name;
  color = kRed;
  if(name == "signal") {
    color = kBlue;
    title = "Signal";
  } else if(name == "mumem") {
    color = kBlue;
    title = "#mu^{-}#rightarrowe^{-}";
  } else if(name == "mumep") {
    color = kBlue;
    title = "#mu^{-}#rightarrowe^{+}";
  } else if(name.BeginsWith("dio")) {
    title = "DIO";
    color = kRed-7;
  } else if(name == "ipa_dio") {
    title = "IPA DIO";
    color = kRed-3;
  } else if(name == "cosmic") {
    title = "Cosmic ray";
    color = kOrange;
  } else if(name == "pbar") {
    title = "Antiproton";
    color = kGreen-6;
  } else if(name == "rpc_ext") {
    title = "RPC (external)";
    color = kMagenta-10;
  } else if(name == "rpc_int") {
    title = "RPC (internal)";
    color = kMagenta+1;
  } else if(name == "rmc_ext") {
    title = "RMC (external)";
    color = kAtlantic+2;
  } else if(name == "rmc_int") {
    title = "RMC (internal)";
    color = kAtlantic;
  } else if(name == "pileup") {
    title = "Pileup";
    color = kGreen-6;
  }
}

#endif
