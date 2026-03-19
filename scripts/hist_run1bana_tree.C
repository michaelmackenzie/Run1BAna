#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TString.h"
#include "TCanvas.h"
#include "TDirectory.h"

#include <string>
#include <vector>
#include <iostream>

//--------------------------------------------------------------------------------------
// Histogram struct - one set per selection
//--------------------------------------------------------------------------------------
struct Hist_t {
  // Cluster
  TH1F* cluster_energy;
  TH1F* cluster_time;
  TH1F* cluster_radius;
  TH1F* cluster_ncr;
  TH1F* cluster_disk;
  TH1F* cluster_e_per_crystal;
  TH1F* cluster_frac_1;
  TH1F* cluster_frac_2;
  TH1F* cluster_second_moment;
  TH1F* cluster_e1;
  TH1F* cluster_e2;
  TH1F* cluster_e9;
  TH1F* cluster_e25;
  TH1F* cluster_t_var;

  // Line-cluster matching
  TH1F* line_dt;
  TH1F* line_dr;
  TH1F* time_cluster_dt;
  TH1F* time_cluster_dr;
  TH1F* ntcl_hits;
  TH1F* photon_id;

  // Line parameters
  TH1F* line_chi2;
  TH1F* line_nhits;
  TH1F* line_nplanes;
  TH1F* line_nstereo;
  TH1F* line_d0;
  TH1F* line_tdip;
  TH1F* line_cos;
  TH1F* line_z0;
  TH1F* line_t0;
  TH1F* line_phi0;

  // Cosmic seed parameters
  TH1F* cosmic_seed_chi2;
  TH1F* cosmic_seed_nhits;
  TH1F* cosmic_seed_d0;
  TH1F* cosmic_seed_tdip;
  TH1F* cosmic_seed_cos;
  TH1F* cosmic_seed_z0;
  TH1F* cosmic_seed_t0;
  TH1F* cosmic_seed_phi0;
  TH1F* cosmic_seed_A0;
  TH1F* cosmic_seed_A1;
  TH1F* cosmic_seed_B0;
  TH1F* cosmic_seed_B1;

  // Time cluster parameters
  TH1F* time_cluster_nhits;
  TH1F* time_cluster_nstraw_hits;
  TH1F* time_cluster_nhigh_z_hits;
  TH1F* time_cluster_t0;
  TH1F* time_cluster_t0err;
  TH1F* time_cluster_z0;
  TH1F* time_cluster_phi0;

  // MC truth
  TH1F* mc_cluster_energy;
  TH1F* mc_cluster_time;
  TH1F* sim_1_edep;
  TH1F* sim_1_time;
  TH1I* sim_1_nhits;
  TH1I* sim_1_type;
  TH1F* sim_2_edep;
  TH1F* sim_2_time;
  TH1I* sim_2_nhits;
  TH1I* sim_2_type;
  TH1I* sim_1_2_nhits;
  TH1F* event_weight;
  TH1F* gen_energy;

  // Directory for this histogram set
  TDirectory* dir = nullptr;
};

//--------------------------------------------------------------------------------------
// Branch variables (mirrors Tree_t in Structs.hh)
//--------------------------------------------------------------------------------------
struct TreeBranches {
  int   event;
  int   subrun;
  int   run;
  float cluster_energy;
  float cluster_time;
  float cluster_radius;
  float cluster_ncr;
  float cluster_disk;
  float cluster_e_per_crystal;
  float cluster_frac_1;
  float cluster_frac_2;
  float cluster_second_moment;
  float cluster_e1;
  float cluster_e2;
  float cluster_e9;
  float cluster_e25;
  float cluster_t_var;
  float line_dt;
  float line_dr;
  float time_cluster_dt;
  float time_cluster_dr;
  float ntcl_hits;
  float photon_id;

  float line_chi2;
  float line_nhits;
  float line_nplanes;
  float line_nstereo;
  float line_d0;
  float line_tdip;
  float line_cos;
  float line_z0;
  float line_t0;
  float line_phi0;

  float cosmic_seed_chi2;
  float cosmic_seed_nhits;
  float cosmic_seed_d0;
  float cosmic_seed_tdip;
  float cosmic_seed_cos;
  float cosmic_seed_z0;
  float cosmic_seed_t0;
  float cosmic_seed_phi0;
  float cosmic_seed_A0;
  float cosmic_seed_A1;
  float cosmic_seed_B0;
  float cosmic_seed_B1;

  float time_cluster_nhits;
  float time_cluster_nstraw_hits;
  float time_cluster_nhigh_z_hits;
  float time_cluster_t0;
  float time_cluster_t0err;
  float time_cluster_z0;
  float time_cluster_phi0;

  float mc_cluster_energy;
  float mc_cluster_time;
  float sim_1_edep;
  float sim_1_time;
  int   sim_1_nhits;
  int   sim_1_type;
  float sim_2_edep;
  float sim_2_time;
  int   sim_2_nhits;
  int   sim_2_type;
  float event_weight;
  float gen_energy;
};

//--------------------------------------------------------------------------------------
// Global histogram array
//--------------------------------------------------------------------------------------
constexpr int kMaxHists = 1000;
Hist_t* hist_[kMaxHists] = {nullptr};

//--------------------------------------------------------------------------------------
double getNSampled(TChain* chain) {
  double total = 0.;
  const TObjArray* files = chain->GetListOfFiles();
  for(int i = 0; i < files->GetEntries(); ++i) {
    const char* fname = files->At(i)->GetTitle();
    TFile* f = TFile::Open(fname, "READ");
    if(!f || f->IsZombie()) {
      std::cerr << "Warning: could not open file " << fname << " for normalization." << std::endl;
      if(f) delete f;
      continue;
    }
    TH1* h = dynamic_cast<TH1*>(f->Get("Run1BAna/data/norm"));
    if(h) total += h->GetEntries();
    else   std::cerr << "Warning: could not retrieve norm histogram from " << fname << std::endl;
    f->Close();
    delete f;
  }
  return total;
}

//--------------------------------------------------------------------------------------
void bookHistograms(const int index, const char* title, TDirectory* outDir) {
  if(index < 0 || index >= kMaxHists) {
    std::cerr << "bookHistograms: index " << index << " out of range!" << std::endl;
    return;
  }
  if(!outDir) {
    std::cerr << "bookHistograms: null output directory" << std::endl;
    return;
  }

  hist_[index] = new Hist_t();
  auto* H = hist_[index];

  const TString setDirName = TString::Format("hist_%d", index);
  H->dir = outDir->mkdir(setDirName, title);
  if(!H->dir) {
    std::cerr << "bookHistograms: failed to create directory " << setDirName << std::endl;
    return;
  }
  H->dir->cd();
  H->dir->SetTitle(title);

  // Cluster
  H->cluster_energy        = new TH1F("cluster_energy"       , "Cluster energy;Energy (MeV);"          , 300,   0.,  300.);
  H->cluster_time          = new TH1F("cluster_time"         , "Cluster time;Time (ns);"               , 200,   0., 2000.);
  H->cluster_radius        = new TH1F("cluster_radius"       , "Cluster radius;R (mm);"                , 100,   0.,  700.);
  H->cluster_ncr           = new TH1F("cluster_ncr"          , "N(crystals);N;"                        ,  20,   0.,   20.);
  H->cluster_disk          = new TH1F("cluster_disk"         , "Disk ID;Disk ID;"                      ,   4,  -1.,    3.);
  H->cluster_e_per_crystal = new TH1F("cluster_e_per_crystal", "Energy/crystal;E/N_{cr} (MeV);"        , 300,   0.,  300.);
  H->cluster_frac_1        = new TH1F("cluster_frac_1"       , "E_{1}/E_{total};Frac;"                 , 101,   0.,  1.01);
  H->cluster_frac_2        = new TH1F("cluster_frac_2"       , "E_{1+2}/E_{total};Frac;"               , 101,   0.,  1.01);
  H->cluster_second_moment = new TH1F("cluster_second_moment", "Second moment;Second moment;"          , 200,   0., 1.e6);
  H->cluster_e1            = new TH1F("cluster_e1"           , "E1;E1 (MeV);"                          , 300,   0.,  300.);
  H->cluster_e2            = new TH1F("cluster_e2"           , "E2;E2 (MeV);"                          , 300,   0.,  300.);
  H->cluster_e9            = new TH1F("cluster_e9"           , "E9;E9 (MeV);"                          , 300,   0.,  300.);
  H->cluster_e25           = new TH1F("cluster_e25"          , "E25;E25 (MeV);"                        , 300,   0.,  300.);
  H->cluster_t_var         = new TH1F("cluster_t_var"        , "Time variance;#sigma_{t}^{2} (ns^{2});", 200,   0.,   10.);

  // Line-cluster matching
  H->line_dt               = new TH1F("line_dt"              , "Line-cluster #Delta t;#Delta t (ns);"  , 200,-200.,  200.);
  H->line_dr               = new TH1F("line_dr"              , "Line-cluster #Delta r;#Delta r (mm);"  , 150,   0.,  500.);
  H->time_cluster_dt       = new TH1F("time_cluster_dt"      , "TCl-cluster #Delta t;#Delta t (ns);"   , 200,-200.,  200.);
  H->time_cluster_dr       = new TH1F("time_cluster_dr"      , "TCl-cluster #Delta r;#Delta r (mm);"   , 150,   0.,  500.);
  H->ntcl_hits             = new TH1F("ntcl_hits"            , "N(TCl hits);N;"                        , 200,   0.,  200.);
  H->photon_id             = new TH1F("photon_id"            , "Photon ID MVA;Score;"                  , 100,  -1.,    1.);

  // Line parameters
  H->line_chi2             = new TH1F("line_chi2"            , "Line #chi^{2}/DOF;#chi^{2}/DOF;"       , 100,   0.,   10.);
  H->line_nhits            = new TH1F("line_nhits"           , "Line N(hits);N;"                       , 100,   0.,  100.);
  H->line_nplanes          = new TH1F("line_nplanes"         , "Line N(planes);N;"                     ,  50,   0.,   50.);
  H->line_nstereo          = new TH1F("line_nstereo"         , "Line N(stereo);N;"                     ,  15,   0.,   15.);
  H->line_d0               = new TH1F("line_d0"              , "Line d_{0};d_{0} (mm);"                , 200,-400.,  400.);
  H->line_tdip             = new TH1F("line_tdip"            , "Line tan(dip);tan(dip);"               , 100, -10.,   10.);
  H->line_cos              = new TH1F("line_cos"             , "Line cos(#theta);cos(#theta);"         , 100,   0.,    1.);
  H->line_z0               = new TH1F("line_z0"              , "Line z_{0};z_{0} (mm);"                , 100,-5000.,5000.);
  H->line_t0               = new TH1F("line_t0"              , "Line t_{0};t_{0} (ns);"                , 200,   0., 2000.);
  H->line_phi0             = new TH1F("line_phi0"            , "Line #phi_{0};#phi_{0} (rad);"         , 100,-3.15,  3.15);

  // Cosmic seed parameters
  H->cosmic_seed_chi2      = new TH1F("cosmic_seed_chi2"     , "Cosmic seed #chi^{2};#chi^{2};"        , 100,   0.,   10.);
  H->cosmic_seed_nhits     = new TH1F("cosmic_seed_nhits"    , "Cosmic seed N(hits);N;"                , 100,   0.,  100.);
  H->cosmic_seed_d0        = new TH1F("cosmic_seed_d0"       , "Cosmic seed d_{0};d_{0} (mm);"         , 200,-400.,  400.);
  H->cosmic_seed_tdip      = new TH1F("cosmic_seed_tdip"     , "Cosmic seed tan(dip);tan(dip);"        , 100, -10.,   10.);
  H->cosmic_seed_cos       = new TH1F("cosmic_seed_cos"      , "Cosmic seed cos(#theta);cos(#theta);"  , 100,   0.,    1.);
  H->cosmic_seed_z0        = new TH1F("cosmic_seed_z0"       , "Cosmic seed z_{0};z_{0} (mm);"         , 100,-5000.,5000.);
  H->cosmic_seed_t0        = new TH1F("cosmic_seed_t0"       , "Cosmic seed t_{0};t_{0} (ns);"         , 200,   0., 2000.);
  H->cosmic_seed_phi0      = new TH1F("cosmic_seed_phi0"     , "Cosmic seed #phi_{0};#phi_{0} (rad);"  , 100,-3.15,  3.15);
  H->cosmic_seed_A0        = new TH1F("cosmic_seed_A0"       , "Cosmic seed A0;A0 (mm);"               , 200,-4000.,4000.);
  H->cosmic_seed_A1        = new TH1F("cosmic_seed_A1"       , "Cosmic seed A1;A1 (mm);"               , 200,-4000.,4000.);
  H->cosmic_seed_B0        = new TH1F("cosmic_seed_B0"       , "Cosmic seed B0;B0 (mm);"               , 200,-4000.,4000.);
  H->cosmic_seed_B1        = new TH1F("cosmic_seed_B1"       , "Cosmic seed B1;B1 (mm);"               , 200,-4000.,4000.);

  // Time cluster parameters
  H->time_cluster_nhits         = new TH1F("time_cluster_nhits"        , "TCl N(hits);N;"              , 200,   0.,  200.);
  H->time_cluster_nstraw_hits   = new TH1F("time_cluster_nstraw_hits"  , "TCl N(straw hits);N;"        , 200,   0.,  200.);
  H->time_cluster_nhigh_z_hits  = new TH1F("time_cluster_nhigh_z_hits" , "TCl N(high-z hits);N;"       ,  50,   0.,   50.);
  H->time_cluster_t0            = new TH1F("time_cluster_t0"           , "TCl t_{0};t_{0} (ns);"       , 200,   0., 2000.);
  H->time_cluster_t0err         = new TH1F("time_cluster_t0err"        , "TCl t_{0} err;t_{0} err (ns);", 100,   0.,   10.);
  H->time_cluster_z0            = new TH1F("time_cluster_z0"           , "TCl z_{0};z_{0} (mm);"       , 100,-5000.,5000.);
  H->time_cluster_phi0          = new TH1F("time_cluster_phi0"         , "TCl #phi_{0};#phi_{0};"       , 100,-3.15,  3.15);

  // MC truth
  H->mc_cluster_energy     = new TH1F("mc_cluster_energy"    , "MC cluster energy;E (MeV);"            , 300,   0.,  300.);
  H->mc_cluster_time       = new TH1F("mc_cluster_time"      , "MC cluster time;t (ns);"               , 200,   0., 2000.);
  H->sim_1_edep            = new TH1F("sim_1_edep"           , "Sim 1 E dep;E (MeV);"                  , 300,   0.,  300.);
  H->sim_1_time            = new TH1F("sim_1_time"           , "Sim 1 time;t (ns);"                    , 200,   0., 2000.);
  H->sim_1_nhits           = new TH1I("sim_1_nhits"          , "Sim 1 N(tracker hits);N;"              , 100,   0,   100);
  H->sim_1_type            = new TH1I("sim_1_type"           , "Sim 1 type;Type;"                      ,  10,  -1,     9);
  H->sim_2_edep            = new TH1F("sim_2_edep"           , "Sim 2 E dep;E (MeV);"                  , 300,   0.,  300.);
  H->sim_2_time            = new TH1F("sim_2_time"           , "Sim 2 time;t (ns);"                    , 200,   0., 2000.);
  H->sim_2_nhits           = new TH1I("sim_2_nhits"          , "Sim 2 N(tracker hits);N;"              , 100,   0,   100);
  H->sim_2_type            = new TH1I("sim_2_type"           , "Sim 2 type;Type;"                      ,  10,  -1,     9);
  H->sim_1_2_nhits         = new TH1I("sim_1_2_nhits"        , "Sim 1-2 N(tracker hits);N;"            , 200,   0,   200);
  H->event_weight          = new TH1F("event_weight"         , "Event weight;Weight;"                  , 100,   0.,    5.);
  H->gen_energy            = new TH1F("gen_energy"           , "Generated energy;E (MeV);"             ,  90,  50.,  140.);
}

//--------------------------------------------------------------------------------------
void fillHistograms(const int index, const TreeBranches& b, double weight = 1.) {
  if(index < 0 || index >= kMaxHists || !hist_[index]) return;
  auto* H = hist_[index];

  const double w = weight;

  // Cluster
  H->cluster_energy       ->Fill(b.cluster_energy,        w);
  H->cluster_time         ->Fill(b.cluster_time,          w);
  H->cluster_radius       ->Fill(b.cluster_radius,        w);
  H->cluster_ncr          ->Fill(b.cluster_ncr,           w);
  H->cluster_disk         ->Fill(b.cluster_disk,          w);
  H->cluster_e_per_crystal->Fill(b.cluster_e_per_crystal, w);
  H->cluster_frac_1       ->Fill(b.cluster_frac_1,        w);
  H->cluster_frac_2       ->Fill(b.cluster_frac_2,        w);
  H->cluster_second_moment->Fill(b.cluster_second_moment, w);
  H->cluster_e1           ->Fill(b.cluster_e1,            w);
  H->cluster_e2           ->Fill(b.cluster_e2,            w);
  H->cluster_e9           ->Fill(b.cluster_e9,            w);
  H->cluster_e25          ->Fill(b.cluster_e25,           w);
  H->cluster_t_var        ->Fill(b.cluster_t_var,         w);

  // Line-cluster matching
  H->line_dt              ->Fill(b.line_dt,               w);
  H->line_dr              ->Fill(b.line_dr,               w);
  H->time_cluster_dt      ->Fill(b.time_cluster_dt,       w);
  H->time_cluster_dr      ->Fill(b.time_cluster_dr,       w);
  H->ntcl_hits            ->Fill(b.ntcl_hits,             w);
  H->photon_id            ->Fill(b.photon_id,             w);

  // Line parameters
  H->line_chi2            ->Fill(b.line_chi2,             w);
  H->line_nhits           ->Fill(b.line_nhits,            w);
  H->line_nplanes         ->Fill(b.line_nplanes,          w);
  H->line_nstereo         ->Fill(b.line_nstereo,          w);
  H->line_d0              ->Fill(b.line_d0,               w);
  H->line_tdip            ->Fill(b.line_tdip,             w);
  H->line_cos             ->Fill(b.line_cos,              w);
  H->line_z0              ->Fill(b.line_z0,               w);
  H->line_t0              ->Fill(b.line_t0,               w);
  H->line_phi0            ->Fill(b.line_phi0,             w);

  // Cosmic seed parameters
  H->cosmic_seed_chi2     ->Fill(b.cosmic_seed_chi2,      w);
  H->cosmic_seed_nhits    ->Fill(b.cosmic_seed_nhits,     w);
  H->cosmic_seed_d0       ->Fill(b.cosmic_seed_d0,        w);
  H->cosmic_seed_tdip     ->Fill(b.cosmic_seed_tdip,      w);
  H->cosmic_seed_cos      ->Fill(b.cosmic_seed_cos,       w);
  H->cosmic_seed_z0       ->Fill(b.cosmic_seed_z0,        w);
  H->cosmic_seed_t0       ->Fill(b.cosmic_seed_t0,        w);
  H->cosmic_seed_phi0     ->Fill(b.cosmic_seed_phi0,      w);
  H->cosmic_seed_A0       ->Fill(b.cosmic_seed_A0,        w);
  H->cosmic_seed_A1       ->Fill(b.cosmic_seed_A1,        w);
  H->cosmic_seed_B0       ->Fill(b.cosmic_seed_B0,        w);
  H->cosmic_seed_B1       ->Fill(b.cosmic_seed_B1,        w);

  // Time cluster parameters
  H->time_cluster_nhits       ->Fill(b.time_cluster_nhits,       w);
  H->time_cluster_nstraw_hits ->Fill(b.time_cluster_nstraw_hits, w);
  H->time_cluster_nhigh_z_hits->Fill(b.time_cluster_nhigh_z_hits,w);
  H->time_cluster_t0          ->Fill(b.time_cluster_t0,          w);
  H->time_cluster_t0err       ->Fill(b.time_cluster_t0err,       w);
  H->time_cluster_z0          ->Fill(b.time_cluster_z0,          w);
  H->time_cluster_phi0        ->Fill(b.time_cluster_phi0,        w);

  // MC truth
  H->mc_cluster_energy    ->Fill(b.mc_cluster_energy,     w);
  H->mc_cluster_time      ->Fill(b.mc_cluster_time,       w);
  H->sim_1_edep           ->Fill(b.sim_1_edep,            w);
  H->sim_1_time           ->Fill(b.sim_1_time,            w);
  H->sim_1_nhits          ->Fill(b.sim_1_nhits,           w);
  H->sim_1_type           ->Fill(b.sim_1_type,            w);
  H->sim_2_edep           ->Fill(b.sim_2_edep,            w);
  H->sim_2_time           ->Fill(b.sim_2_time,            w);
  H->sim_2_nhits          ->Fill(b.sim_2_nhits,           w);
  H->sim_2_type           ->Fill(b.sim_2_type,            w);
  H->sim_1_2_nhits        ->Fill(b.sim_1_nhits + b.sim_2_nhits, w);
  H->event_weight         ->Fill(b.event_weight,          w);
  H->gen_energy           ->Fill(b.gen_energy,            w);
}

//--------------------------------------------------------------------------------------
void setBranchAddresses(TTree* tree, TreeBranches& b) {
  tree->SetBranchAddress("event"                   , &b.event);
  tree->SetBranchAddress("subrun"                  , &b.subrun);
  tree->SetBranchAddress("run"                     , &b.run);
  tree->SetBranchAddress("cluster_energy"          , &b.cluster_energy);
  tree->SetBranchAddress("cluster_time"            , &b.cluster_time);
  tree->SetBranchAddress("cluster_radius"          , &b.cluster_radius);
  tree->SetBranchAddress("cluster_ncr"             , &b.cluster_ncr);
  tree->SetBranchAddress("cluster_disk"            , &b.cluster_disk);
  tree->SetBranchAddress("cluster_e_per_crystal"   , &b.cluster_e_per_crystal);
  tree->SetBranchAddress("cluster_frac_1"          , &b.cluster_frac_1);
  tree->SetBranchAddress("cluster_frac_2"          , &b.cluster_frac_2);
  tree->SetBranchAddress("cluster_second_moment"   , &b.cluster_second_moment);
  tree->SetBranchAddress("cluster_e1"              , &b.cluster_e1);
  tree->SetBranchAddress("cluster_e2"              , &b.cluster_e2);
  tree->SetBranchAddress("cluster_e9"              , &b.cluster_e9);
  tree->SetBranchAddress("cluster_e25"             , &b.cluster_e25);
  tree->SetBranchAddress("cluster_t_var"           , &b.cluster_t_var);
  tree->SetBranchAddress("line_dt"                 , &b.line_dt);
  tree->SetBranchAddress("line_dr"                 , &b.line_dr);
  tree->SetBranchAddress("time_cluster_dt"         , &b.time_cluster_dt);
  tree->SetBranchAddress("time_cluster_dr"         , &b.time_cluster_dr);
  tree->SetBranchAddress("ntcl_hits"               , &b.ntcl_hits);
  tree->SetBranchAddress("photon_id"               , &b.photon_id);
  tree->SetBranchAddress("line_chi2"               , &b.line_chi2);
  tree->SetBranchAddress("line_nhits"              , &b.line_nhits);
  tree->SetBranchAddress("line_nplanes"            , &b.line_nplanes);
  tree->SetBranchAddress("line_nstereo"            , &b.line_nstereo);
  tree->SetBranchAddress("line_d0"                 , &b.line_d0);
  tree->SetBranchAddress("line_tdip"               , &b.line_tdip);
  tree->SetBranchAddress("line_cos"                , &b.line_cos);
  tree->SetBranchAddress("line_z0"                 , &b.line_z0);
  tree->SetBranchAddress("line_t0"                 , &b.line_t0);
  tree->SetBranchAddress("line_phi0"               , &b.line_phi0);
  tree->SetBranchAddress("cosmic_seed_chi2"        , &b.cosmic_seed_chi2);
  tree->SetBranchAddress("cosmic_seed_nhits"       , &b.cosmic_seed_nhits);
  tree->SetBranchAddress("cosmic_seed_d0"          , &b.cosmic_seed_d0);
  tree->SetBranchAddress("cosmic_seed_tdip"        , &b.cosmic_seed_tdip);
  tree->SetBranchAddress("cosmic_seed_cos"         , &b.cosmic_seed_cos);
  tree->SetBranchAddress("cosmic_seed_z0"          , &b.cosmic_seed_z0);
  tree->SetBranchAddress("cosmic_seed_t0"          , &b.cosmic_seed_t0);
  tree->SetBranchAddress("cosmic_seed_phi0"        , &b.cosmic_seed_phi0);
  tree->SetBranchAddress("cosmic_seed_A0"          , &b.cosmic_seed_A0);
  tree->SetBranchAddress("cosmic_seed_A1"          , &b.cosmic_seed_A1);
  tree->SetBranchAddress("cosmic_seed_B0"          , &b.cosmic_seed_B0);
  tree->SetBranchAddress("cosmic_seed_B1"          , &b.cosmic_seed_B1);
  tree->SetBranchAddress("time_cluster_nhits"      , &b.time_cluster_nhits);
  tree->SetBranchAddress("time_cluster_nstraw_hits", &b.time_cluster_nstraw_hits);
  tree->SetBranchAddress("time_cluster_nhigh_z_hits",&b.time_cluster_nhigh_z_hits);
  tree->SetBranchAddress("time_cluster_t0"         , &b.time_cluster_t0);
  tree->SetBranchAddress("time_cluster_t0err"      , &b.time_cluster_t0err);
  tree->SetBranchAddress("time_cluster_z0"         , &b.time_cluster_z0);
  tree->SetBranchAddress("time_cluster_phi0"       , &b.time_cluster_phi0);
  tree->SetBranchAddress("mc_cluster_energy"       , &b.mc_cluster_energy);
  tree->SetBranchAddress("mc_cluster_time"         , &b.mc_cluster_time);
  tree->SetBranchAddress("sim_1_edep"              , &b.sim_1_edep);
  tree->SetBranchAddress("sim_1_time"              , &b.sim_1_time);
  tree->SetBranchAddress("sim_1_nhits"             , &b.sim_1_nhits);
  tree->SetBranchAddress("sim_1_type"              , &b.sim_1_type);
  tree->SetBranchAddress("sim_2_edep"              , &b.sim_2_edep);
  tree->SetBranchAddress("sim_2_time"              , &b.sim_2_time);
  tree->SetBranchAddress("sim_2_nhits"             , &b.sim_2_nhits);
  tree->SetBranchAddress("sim_2_type"              , &b.sim_2_type);
  tree->SetBranchAddress("event_weight"            , &b.event_weight);
  tree->SetBranchAddress("gen_energy"              , &b.gen_energy);
}

//--------------------------------------------------------------------------------------
// Selection functions
//--------------------------------------------------------------------------------------

bool sel_energy(const TreeBranches& b) {
  return b.cluster_energy > 60.;
}

bool sel_energy_time(const TreeBranches& b) {
  return b.cluster_energy > 60. && b.cluster_time > 400. && b.cluster_time < 1650.;
}

bool sel_photon_id(const TreeBranches& b) {
  return sel_energy_time(b) && b.photon_id > 0.8;
}

bool sel_signal_id(const TreeBranches& b) {
  return sel_energy_time(b)
    && b.cluster_ncr  > 1
    && b.cluster_ncr  < 6
                        && b.cluster_frac_1       > 0.60f
    && b.cluster_frac_2       > 0.80f
    && b.cluster_t_var        < 1.0f
    && b.cluster_second_moment< 1.e5f;
  // && b.photon_id            > 0.8f;
}

//--------------------------------------------------------------------------------------
// Main entry point
//--------------------------------------------------------------------------------------
void hist_run1bana_tree(const char* inputFiles = "input.root",  // comma- or space-separated, or a glob
                        const char* outputFile = "output.root",
                        const char* treePath   = "Run1BAna/tree_60/tree") {

  // Build TChain
  TChain* chain = new TChain(treePath);

  // Allow comma-separated list of files or wildcards
  TString fileList(inputFiles);
  TObjArray* tokens = fileList.Tokenize(",");
  for(int i = 0; i < tokens->GetEntries(); ++i) {
    TString fname = dynamic_cast<TObjString*>(tokens->At(i))->GetString().Strip(TString::kBoth);
    const int added = chain->Add(fname);
    std::cout << "Added " << added << " file(s) matching: " << fname << std::endl;
  }
  delete tokens;

  if(chain->GetNtrees() == 0) {
    std::cerr << "No files added to chain, exiting." << std::endl;
    delete chain;
    return;
  }

  // Get normalization by scanning input files before the event loop
  const double nsampled = getNSampled(chain);
  if(nsampled > 0.) {
    std::cout << "Total number of sampled events: " << nsampled << std::endl;
  } else {
    std::cerr << "Warning: could not retrieve number of sampled events from input files." << std::endl;
    delete chain;
    return;
  }

  // Initialise branch addresses
  TreeBranches b{};
  setBranchAddresses(chain, b);

  // Create output file
  TFile* fout = TFile::Open(outputFile, "RECREATE");
  if(!fout || fout->IsZombie()) {
    std::cerr << "Cannot create output file: " << outputFile << std::endl;
    delete chain;
    return;
  }

  // Add normalization to the output
  fout->cd();
  TH1* hnorm = new TH1D("norm", "Normalization;N;", 1, 0., 1.);
  hnorm->SetBinContent(1, nsampled);
  hnorm->Write();

  // Book histogram sets
  bookHistograms(  0, "all"                                   , fout);
  bookHistograms(  1, "no_weights"                            , fout);
  bookHistograms(  2, "photon_id"                             , fout);
  bookHistograms(  3, "signal_id"                             , fout);
  bookHistograms(  4, "id_high_z_hits"                        , fout);
  bookHistograms( 10, "r_500"                                 , fout);
  bookHistograms( 11, "r_550"                                 , fout);
  bookHistograms( 15, "id_r_500"                              , fout);
  bookHistograms( 16, "id_r_550"                              , fout);
  bookHistograms( 17, "id_no_hits"                            , fout);
  bookHistograms( 18, "id_r_500_high_z_hits"                  , fout);
  bookHistograms( 20, "no_calo_mu"                            , fout);
  bookHistograms( 22, "no_calo_mu_photon_id"                  , fout);
  bookHistograms( 23, "no_calo_mu_id"                         , fout);
  bookHistograms( 24, "no_calo_mu_id_high_z_hits"             , fout);
  bookHistograms( 30, "no_calo_mu_r_500"                      , fout);
  bookHistograms( 31, "no_calo_mu_r_550"                      , fout);
  bookHistograms( 32, "n_calo_mu_no_hits"                     , fout);
  bookHistograms( 35, "no_calo_mu_id_r_500"                   , fout);
  bookHistograms( 36, "no_calo_mu_id_r_550"                   , fout);
  bookHistograms( 37, "no_calo_mu_id_no_hits"                 , fout);
  bookHistograms( 38, "no_calo_mu_id_r_500_high_z_hits"       , fout);

  // Sets with offsets
  for(int offset = 0; offset < 3; ++offset) {

    // RMC sets
    bookHistograms( 70 + offset*100, "base"             , fout);
    bookHistograms( 71 + offset*100, "id"               , fout);
    bookHistograms( 72 + offset*100, "id_r_500"         , fout);
    bookHistograms( 73 + offset*100, "id_tcl_hits"      , fout);
    bookHistograms( 74 + offset*100, "id_r_500_tcl_hits", fout);

    // RPC sets
    bookHistograms( 90 + offset*100, "base"             , fout);
    bookHistograms( 91 + offset*100, "id"               , fout);
    bookHistograms( 92 + offset*100, "id_r_500"         , fout);
    bookHistograms( 93 + offset*100, "id_tcl_hits"      , fout);
    bookHistograms( 94 + offset*100, "id_r_500_tcl_hits", fout);

    bookHistograms( 95 + offset*100, "t_500", fout);
    bookHistograms( 96 + offset*100, "id_t_500", fout);
    bookHistograms( 97 + offset*100, "sim_t_500", fout);
  }


  //--------------------------------------------------------------------------------------
  // Event loop
  //--------------------------------------------------------------------------------------
  const Long64_t nEntries = chain->GetEntries();
  std::cout << "Processing " << nEntries << " entries across "
            << chain->GetNtrees() << " file(s)..." << std::endl;

  const bool is_pu = TString(inputFiles).Contains("mnbs");

  for(Long64_t i = 0; i < nEntries; ++i) {
    chain->GetEntry(i);

    if(i % 10000 == 0)
      std::cout << "  Entry " << i << " / " << nEntries << std::endl;

    int offset = 0;
    if(b.sim_1_type == 2 || b.sim_2_type == 2) offset = 200; // calo muon stop
    else if(b.sim_1_edep / b.cluster_energy < 0.75) offset = 100; // pileup or misreconstructed

    // Fill each selection set
    fillHistograms(0, b, b.event_weight);
    fillHistograms(1, b);
    if(sel_photon_id    (b)) fillHistograms(2, b, b.event_weight);
    if(sel_signal_id    (b)) {
      fillHistograms(3, b, b.event_weight);
      if(b.time_cluster_nhigh_z_hits < 3) fillHistograms(4, b, b.event_weight);
    }
    if(sel_energy_time(b)) {
      if(b.cluster_radius > 500.)  fillHistograms(10, b, b.event_weight);
      if(b.cluster_radius > 550.)  fillHistograms(11, b, b.event_weight);
      if(sel_signal_id(b)) {
        if(b.cluster_radius > 500.) {
          fillHistograms(15, b, b.event_weight);
          if(b.time_cluster_nhigh_z_hits < 3) fillHistograms(18, b, b.event_weight);
        }
        if(b.cluster_radius > 550.)  fillHistograms(16, b, b.event_weight);
        if(b.sim_1_nhits + b.sim_2_nhits <= 0) fillHistograms(17, b, b.event_weight);
      }
    }

    // No calorimter muon stops
    if(b.sim_1_type != 2 && b.sim_2_type != 2) {
      fillHistograms(20, b, b.event_weight);
      if(sel_photon_id    (b)) fillHistograms(22, b, b.event_weight);
      if(sel_signal_id    (b)) {
        fillHistograms(23, b, b.event_weight);
        if(b.time_cluster_nhigh_z_hits < 3) fillHistograms(24, b, b.event_weight);
      }
      if(sel_energy_time(b)) {
        if(b.cluster_radius > 500.)  fillHistograms(30, b, b.event_weight);
        if(b.cluster_radius > 550.)  fillHistograms(31, b, b.event_weight);
        if(b.sim_1_nhits + b.sim_2_nhits <= 0) fillHistograms(32, b, b.event_weight);
        if(sel_signal_id(b)) {
          if(b.cluster_radius > 500.)  {
            fillHistograms(35, b, b.event_weight);
            if(b.time_cluster_nhigh_z_hits < 3) fillHistograms(38, b, b.event_weight);
          }
          if(b.cluster_radius > 550.)  fillHistograms(36, b, b.event_weight);
          if(b.sim_1_nhits + b.sim_2_nhits <= 0) fillHistograms(37, b, b.event_weight);
        }
      }
    }

    // Offset selections

    // RMC
    if(b.cluster_time > 600.) {
      fillHistograms(70 + offset, b, b.event_weight);
      if(sel_signal_id(b)) {
        fillHistograms(71 + offset, b, b.event_weight);
        if(b.cluster_radius > 500.) fillHistograms(72 + offset, b, b.event_weight);
        if(b.time_cluster_nhigh_z_hits < 3) fillHistograms(73 + offset, b, b.event_weight);
        if(b.cluster_radius > 500. && b.time_cluster_nhigh_z_hits < 3) fillHistograms(74 + offset, b, b.event_weight);
      }
    }

    // Low time selection
    if(b.cluster_time > 400. && b.cluster_time < 550.) {
      fillHistograms(90 + offset, b, b.event_weight);
      const bool signal_id = (   b.cluster_energy > 60.
                              && b.cluster_ncr  > 2
                              && b.cluster_ncr  < 8
                              // && b.cluster_frac_1       > 0.60f
                              // && b.cluster_frac_2       > 0.80f
                              && b.cluster_t_var        < 1.f
                              && b.cluster_second_moment< 2.e5f
                              );
      if(signal_id) {
        fillHistograms(91 + offset, b, b.event_weight);
        if(b.cluster_radius > 500.) fillHistograms(92 + offset, b, b.event_weight);
        if(b.time_cluster_nhigh_z_hits < 3) fillHistograms(93 + offset, b, b.event_weight);
        if(b.cluster_radius > 500. && b.time_cluster_nhigh_z_hits < 3) fillHistograms(94 + offset, b, b.event_weight);
      }
      if(b.cluster_time > 500.) {
        fillHistograms(95 + offset, b, b.event_weight);
        if(signal_id) fillHistograms(96 + offset, b, b.event_weight);
      }
      if(b.sim_1_time > 500.) fillHistograms(97 + offset, b, b.event_weight);
    }

    // if(offset == 0 && is_pu) { // report relevant event IDs
    //   std::cout << "PU DIO: " << b.run << " " << b.subrun << " " << b.event << std::endl;
    // }
  }

  //--------------------------------------------------------------------------------------
  // Save output
  //--------------------------------------------------------------------------------------
  fout->Write();
  fout->Close();
  delete chain;

  std::cout << "Done. Output written to " << outputFile << std::endl;
}
