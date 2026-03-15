#include "TFile.h"
#include "TTree.h"
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
  TH1F* event_weight;
  TH1F* gen_energy;

  // Directory for this histogram set
  TDirectory* dir = nullptr;
};

//--------------------------------------------------------------------------------------
// Branch variables (mirrors Tree_t in Structs.hh)
//--------------------------------------------------------------------------------------
struct TreeBranches {
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
constexpr int kMaxHists = 20;
Hist_t* hist_[kMaxHists] = {nullptr};

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
  H->dir = outDir->mkdir(setDirName);
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
  H->event_weight          = new TH1F("event_weight"         , "Event weight;Weight;"                  , 100,   0.,    5.);
  H->gen_energy            = new TH1F("gen_energy"           , "Generated energy;E (MeV);"             ,  60,  50.,  110.);
}

//--------------------------------------------------------------------------------------
void fillHistograms(const int index, const TreeBranches& b, double weight = 1.) {
  if(index < 0 || index >= kMaxHists || !hist_[index]) return;
  auto* H = hist_[index];

  const double w = b.event_weight * weight;

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
  H->event_weight         ->Fill(b.event_weight,          w);
  H->gen_energy           ->Fill(b.gen_energy,            w);
}

//--------------------------------------------------------------------------------------
void setBranchAddresses(TTree* tree, TreeBranches& b) {
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
bool sel_all(const TreeBranches& /*b*/) {
  return true;
}

bool sel_70MeV(const TreeBranches& b) {
  return b.cluster_energy > 70.;
}

bool sel_70MeV_in_time(const TreeBranches& b) {
  return b.cluster_energy > 70. && b.cluster_time > 600. && b.cluster_time < 1650.;
}

bool sel_photon_id(const TreeBranches& b) {
  return sel_70MeV_in_time(b) && b.photon_id > 0.8;
}

bool sel_signal_id(const TreeBranches& b) {
  return sel_70MeV_in_time(b)
      && b.cluster_ncr  > 1
      && b.cluster_ncr  < 6
      && b.cluster_frac_1       > 0.60f
      && b.cluster_frac_2       > 0.80f
      && b.cluster_t_var        < 1.0f
      && b.cluster_second_moment< 1.e5f
      && b.photon_id            > 0.8f;
}

//--------------------------------------------------------------------------------------
// Main entry point
//--------------------------------------------------------------------------------------
void analyze_tree(const char* inputFile  = "input.root",
                  const char* treePath   = "tree_60/tree",
                  const char* outputFile = "output.root")
{
  // Create output file first so booking can create subdirectories immediately
  TFile* fout = TFile::Open(outputFile, "RECREATE");
  if(!fout || fout->IsZombie()) {
    std::cerr << "Cannot create output file: " << outputFile << std::endl;
    return;
  }

  // Book histogram sets (each into hist_<index>)
  bookHistograms(0, "all",          fout);
  bookHistograms(1, "70MeV",        fout);
  bookHistograms(2, "70MeV_in_time",fout);
  bookHistograms(3, "photon_id",    fout);
  bookHistograms(4, "signal_id",    fout);

  // Open input file and tree
  TFile* fin = TFile::Open(inputFile, "READ");
  if(!fin || fin->IsZombie()) {
    std::cerr << "Cannot open input file: " << inputFile << std::endl;
    fout->Close();
    return;
  }

  // Initialise branch addresses
  TreeBranches b{};
  setBranchAddresses(tree, b);

  const Long64_t nEntries = tree->GetEntries();
  std::cout << "Processing " << nEntries << " entries..." << std::endl;

  //--------------------------------------------------------------------------------------
  // Event loop
  //--------------------------------------------------------------------------------------
  for(Long64_t i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);

    if(i % 10000 == 0)
      std::cout << "  Entry " << i << " / " << nEntries << std::endl;

    // Fill each selection set
    if(sel_all          (b)) fillHistograms(0, b);
    if(sel_70MeV        (b)) fillHistograms(1, b);
    if(sel_70MeV_in_time(b)) fillHistograms(2, b);
    if(sel_photon_id    (b)) fillHistograms(3, b);
    if(sel_signal_id    (b)) fillHistograms(4, b);
  }

  //--------------------------------------------------------------------------------------
  // Save output
  //--------------------------------------------------------------------------------------
  writeHistograms(fout);
  fout->Write();
  fout->Close();
  fin->Close();

  std::cout << "Done. Output written to " << outputFile << std::endl;
}