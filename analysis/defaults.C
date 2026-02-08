#ifndef __RUN1BANA_ANALYSIS_DEFAULTS__
#define __RUN1BANA_ANALYSIS_DEFAULTS__

const char*  hist_path_   = "/exp/mu2e/data/users/mmackenz/run1b/histograms/";
TString      var_         = "obs"; //observable name
const double bin_width_   =  0.2;  //expected bin width, rebin to achieve if possible
const bool   hist_pdfs_   = true;  //use functions or histograms in the model
const bool   include_sys_ = true;  //evaluate systematics

#endif
