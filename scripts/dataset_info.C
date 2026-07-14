// Define input dataset information

#ifndef RUN1BANA_DATASET_INFO_C
#define RUN1BANA_DATASET_INFO_C

#include "Run1BAna/analysis/physics.C"

struct Dataset_t {
  TString name;
  Double_t nGen;  // N(generated events) or livetime for cosmic datasets
  Double_t nDigi; // N(digitized events) in the parent dataset
  Double_t xsec;  // Rate per POT for beam processes
  Double_t emin;  // Generation range for flat datasets
  Double_t emax;
  TString parent; // Input digi dataset name

  Dataset_t(TString name_ = "", Double_t nGen_ = 0, Double_t nDigi_ = 0, Double_t xsec_ = 0, TString parent_ = "")
    : name(name_), nGen(nGen_), nDigi(nDigi_), xsec(xsec_), emin(0.), emax(1.), parent(parent_) {}
  Dataset_t(TString name_, Double_t nGen_, Double_t nDigi_, Double_t xsec_, Double_t emin_, Double_t emax_, TString parent_ = "")
    : name(name_), nGen(nGen_), nDigi(nDigi_), xsec(xsec_), emin(emin_), emax(emax_), parent(parent_) {}

  bool isCosmic() const { return name.BeginsWith("csms"); }
  double norm(double npot, double livetime) const {
    if (isCosmic()) {
      return livetime / nGen;
    } else {
      return (npot * xsec * (emax - emin)) / nGen;
    }
  }
  TString fileName() const {
    TString histName = "Run1BAna." + name + ".hist";
    return histName;
  }

};

// Retrieve known datasets
map<TString, Dataset_t> getDatasets(TString version = "v40") {
  map<TString, Dataset_t> datasets;

  // Rates without configuration-specific stopping rates
  const double rate_ce        = muon_capture_fraction_;
  const double rate_neut_calo = muon_capture_fraction_*1.2*5.05038e-06; // N(neutrons per capture)*(KE > 60)*(cz restriction)
  const double rate_prot_calo = muon_capture_fraction_*0.05*1.10588e-05; // N(protons per capture)*(KE > 60)*(cz restriction)
  const double rate_dio       = (1. - muon_capture_fraction_);
  const double rate_dio_50    = rate_dio*dio_frac_50_;
  const double rate_dio_80    = rate_dio*dio_frac_80_;
  const double rate_dio_90    = rate_dio*dio_frac_90_;
  const double rate_dio_95    = rate_dio*dio_frac_95_;
  const double rate_dio_0_60  = rate_dio*dio_frac_0_60_;
  const double rate_dio_60_80 = rate_dio*dio_frac_60_80_;
  const double rate_dio_80_90 = rate_dio*dio_frac_80_90_;
  const double rate_rpc       = pion_stop_rate_*pion_survive_frac_*rpc_br_; // FIXME
  const double rate_rpc_int   = rate_rpc*rpc_int_br_; // FIXME
  const double rate_rmc       = muon_capture_fraction_*br_rmc_; // inputs are normalized above 57 MeV
  const double rate_rmc_conv  = rate_rmc*rmc_conv_;
  const double rate_rmc_int   = rate_rmc*rmc_int_br_;

  if(version == "v40") {
    const double muons_per_pot = 5.066e-04;
    nmuons_per_pot_ = muons_per_pot;
    datasets.emplace("mnbs", Dataset_t("mnbs0b1s51r0003",   99995000, 99995000, 1.0                         ,            "dig.mu2e.NoPrimaryMix1BB.Run1Ban_best_v1_4-000.art"));
    datasets.emplace("cele", Dataset_t("cele0b1s51r0003", 1999000000,  1326786, rate_ce *muons_per_pot      ,            "dig.mu2e.CeEndpointMix1BB.Run1Ban_best_v1_4-000.art"));
    datasets.emplace("fgam", Dataset_t("fgam0b1s51r0003", 1999000000,  1039674, rate_rmc*muons_per_pot      , 50., 110., "dig.mu2e.FlatGammaMix1BB.Run1Ban_best_v1_4-000.art"));
    datasets.emplace("csms", Dataset_t("csms0b1s51r0003",    18428.5,  2351533, 1.0                         ,            "dig.mu2e.CosmicCRYAllMix1BB.Run1Ban_best_v1_4-000.art"));
    datasets.emplace("dio0", Dataset_t("dio00b1s51r0003",         1.,       1., 1.0                         ,            "dig.mu2e.DIOMix1BB.Run1Ban_best_v1_4-000.art"));
    datasets.emplace("neut", Dataset_t("neut0b1s51r0003",    3450000,    16704, rate_neut_calo*muons_per_pot,            "dig.mu2e.neut0b0s41r0000.Run1Ban_best_v1_4-000.art"));
    datasets.emplace("prot", Dataset_t("prot0b1s51r0003",     100000,      633, rate_prot_calo*muons_per_pot,            "dig.mu2e.prot0b0s41r0000.Run1Ban_best_v1_4-000.art"));
  } else {
    std::cerr << "Unknown dataset version: " << version << std::endl;
  }
  return datasets;
}

// Validate dataset normalization information
void validateDatasets(TString version = "v40") {
  auto datasets = getDatasets(version);
  for (const auto& [name, ds] : datasets) {

    // Test the digi count info
    TString command = "./Run1BAna/scripts/samCountEvents.sh " + ds.parent + " || echo 0";
    long nDigi = 0;
    FILE* pipe = popen(command.Data(), "r");
    if (pipe) {
      fscanf(pipe, "%ld", &nDigi);
      pclose(pipe);
    }
    if (nDigi != ds.nDigi) {
      std::cerr << "Warning: Dataset " << ds.name << " has nDigi = " << nDigi
                << ", expected " << ds.nDigi << std::endl;
    }

    // Test the (non-cosmic) N(gen events) info
    command = "./Run1BAna/scripts/samCountGenEvents.sh " + ds.parent + " | tail -n 1 | awk -F'=' '{print $NF}' || echo 0";
    long nGen = 0;
    pipe = popen(command.Data(), "r");
    if (pipe) {
      fscanf(pipe, "%ld", &nGen);
      pclose(pipe);
    }
    if(!ds.isCosmic() && nGen != ds.nGen) {
      std::cerr << "Warning: Dataset " << ds.name << " has nGen = " << nGen
                << ", expected " << ds.nGen << std::endl;
    }
    printf("Dataset %s: nGen = %12ld, nDigi = %12ld, xsec = %12.3g, parent = %s\n",
           ds.name.Data(), nGen, nDigi, ds.xsec, ds.parent.Data());
  }
}

#endif // RUN1BANA_DATASET_INFO_C
