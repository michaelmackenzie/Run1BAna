// Useful structs

// Offline
#include "Offline/CaloCluster/inc/ClusterUtils.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

// ROOT
#include "TH1.h"
#include "TH2.h"

using namespace mu2e;

namespace Run1BAnaStructs {

  //--------------------------------------------------------------------------------------
  // Histograms
  //--------------------------------------------------------------------------------------

  // Per event
  struct EventHist_t {
    TH1* npot;
    TH1* n_mc_digis;
    TH1* ncombo_hits;
    TH1* ncalo_hits;
    TH1* nclusters;
    TH1* ngood_clusters;
    TH1* nlines;
    TH1* ngood_lines;
    TH1* ncosmic_seeds;
    TH1* ngood_cosmic_seeds;
    TH1* ntime_clusters;
    TH1* ngood_time_clusters;
    TH1* trig_bits;
    TH1* trig_paths;
    TH2* sim_dr_dt;
    TH2* hit_x_y;
  };

      // Per cluster
    struct ClusterHist_t {
      TH1* energy;
      TH1* mc_edep;
      TH1* time;
      TH1* radius;
      TH1* ncr;
      TH1* disk;
      TH1* energy_per_crystal;
      TH1* frac_first_crystal;
      TH1* frac_first_two_crystals;
      TH1* second_moment;
      TH1* e1;
      TH1* e2;
      TH1* e9;
      TH1* e25;
      TH1* t_var;
      TH1* photon_id;

      TH1* line_dt;
      TH1* line_dr;
      TH1* nmatched_lines;
      TH1* nfit_matched_lines;

      TH1* time_cluster_dt;
      TH1* time_cluster_dr;
      TH1* nmatched_time_clusters;
      TH1* nfit_matched_time_clusters;

      TH1* esum;
      TH1* dt;
      TH1* dr;
      TH2* dr_vs_dt;

      TH1* pdg;
      TH1* energy_sim;
      TH1* energy_ratio;
      TH1* sim_1_nhits;
      TH1* sim_1_type;
      TH1* pdg2;
      TH1* energy_sim2;
      TH1* energy_ratio2;
      TH1* sim_2_nhits;
      TH1* sim_2_type;
      TH1* sim_1_2_nhits;
      TH1* sim_dt;
      TH1* sim_dr;
      TH2* sim_dt_dr;
      TH1* MC_energy_diff;
      TH1* MC_time_diff;
      TH2* energy_time_diff2d;
      TH2* energy_vs_gen_energy;
    };

    // Per line (KalSeed)
    struct LineHist_t {
      TH1* chi2;
      TH1* nhits;
      TH1* nplanes;
      TH1* nstereo;
      TH1* d0;
      TH1* tdip;
      TH1* cos;
      TH1* z0;
      TH1* t0;
      TH1* phi0;

      TH1* cl_energy; // associated calo cluster info
      TH1* cl_dt;
      TH1* cl_dr;
    };

    // Per cosmic seed
    struct CosmicSeedHist_t {
      TH1* chi2;
      TH1* nhits;
      TH1* d0;
      TH1* tdip;
      TH1* cos;
      TH1* z0;
      TH1* t0;
      TH1* phi0;

      TH1* A0;
      TH1* A1;
      TH1* B0;
      TH1* B1;

      TH1* cl_energy; // associated calo cluster info
      TH1* cl_time;
      TH1* cl_disk;
      TH1* cl_dt;
      TH1* cl_dr;
    };

    // Per time cluster
    struct TimeClusterHist_t {
      TH1* nhits;
      TH1* nstraw_hits;
      TH1* nhigh_z_hits;
      TH1* t0;
      TH1* t0err;
      TH1* z0;
      TH1* phi0; // phi in the tracker system
      TH1* cl_energy;
      TH1* cl_time;
      TH1* cl_disk;
      TH1* cl_dt;
      TH1* cl_dr;
      TH1* n_primary_hits;
      TH1* n_other_hits;
      TH1* purity;
      TH1* efficiency;
    };

    // Per sim
    struct SimHist_t {
      TH1* pdg;
      TH1* type;
      TH1* parent_type;
      TH1* origin_type;
      TH1* nhits;
      TH1* energy_start;
      TH1* time_start;
      TH1* mom_cz_start;
      TH1* time_end;
      TH1* edep;
      TH2* start_x_y;
      TH2* energy_vs_trig_path;
    };

    // Output tree branches
    struct Tree_t {

      // Event info
      int   event;
      int   subrun;
      int   run;
      float event_weight;

      // Cluster info
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

      // Line info
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

      // Cosmic seed info
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

      // Time cluster info
      float time_cluster_nhits;
      float time_cluster_nstraw_hits;
      float time_cluster_nhigh_z_hits;
      float time_cluster_t0;
      float time_cluster_t0err;
      float time_cluster_z0;
      float time_cluster_phi0;

      // MC truth info
      float mc_cluster_energy;
      float mc_cluster_time;
      float sim_1_edep;
      float sim_1_time;
      int   sim_1_nhits;
      int   sim_1_type;
      int   sim_1_pdg;
      int   sim_1_main_crystal;
      float sim_1_main_crystal_energy;
      float sim_2_edep;
      float sim_2_time;
      int   sim_2_nhits;
      int   sim_2_type;
      int   sim_2_pdg;
      int   sim_2_main_crystal;
      float sim_2_main_crystal_energy;
      float gen_energy;
      float npot;

      Tree_t() {
        init();
      }

      void init() {
        event_weight = 1.f;
        event = 0;
        subrun = 0;
        run = 0;

        cluster_energy = 0.f;
        cluster_time = 0.f;
        cluster_radius = 0.f;
        cluster_ncr = 0.f;
        cluster_disk = -1.f;
        cluster_e_per_crystal = 0.f;
        cluster_frac_1 = 0.f;
        cluster_frac_2 = 0.f;
        cluster_second_moment = 0.f;
        cluster_e1 = 0.f;
        cluster_e2 = 0.f;
        cluster_e9 = 0.f;
        cluster_e25 = 0.f;
        cluster_t_var = 0.f;
        line_dt = 0.f;
        line_dr = 0.f;
        time_cluster_dt = 0.f;
        time_cluster_dr = 0.f;
        ntcl_hits = 0.f;
        photon_id = 0.f;

        line_chi2 = 0.f;
        line_nhits = 0.f;
        line_nplanes = 0.f;
        line_nstereo = 0.f;
        line_d0 = 0.f;
        line_tdip = 0.f;
        line_cos = 0.f;
        line_z0 = 0.f;
        line_t0 = 0.f;
        line_phi0 = 0.f;

        cosmic_seed_chi2 = 0.f;
        cosmic_seed_nhits = 0.f;
        cosmic_seed_d0 = 0.f;
        cosmic_seed_tdip = 0.f;
        cosmic_seed_cos = 0.f;
        cosmic_seed_z0 = 0.f;
        cosmic_seed_t0 = 0.f;
        cosmic_seed_phi0 = 0.f;
        cosmic_seed_A0 = 0.f;
        cosmic_seed_A1 = 0.f;
        cosmic_seed_B0 = 0.f;
        cosmic_seed_B1 = 0.f;

        time_cluster_nhits = 0.f;
        time_cluster_nstraw_hits = 0.f;
        time_cluster_nhigh_z_hits = 0.f;
        time_cluster_t0 = 0.f;
        time_cluster_t0err = 0.f;
        time_cluster_z0 = 0.f;
        time_cluster_phi0 = 0.f;

        mc_cluster_energy = 0.f;
        mc_cluster_time = 0.f;
        sim_1_edep = 0.f;
        sim_1_time = 0.f;
        sim_1_nhits = 0;
        sim_1_type = -1;
        sim_1_pdg = 0;
        sim_1_main_crystal = -1;
        sim_1_main_crystal_energy = 0.f;
        sim_2_edep = 0.f;
        sim_2_time = 0.f;
        sim_2_nhits = 0;
        sim_2_type = -1;
        sim_2_pdg = 0;
        sim_2_main_crystal = -1;
        sim_2_main_crystal_energy = 0.f;
        gen_energy = 0.f;
        npot = 0.f;
      }
    };


    //--------------------------------------------------------------------------------------
    // Internal data structures
    //--------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------
    struct EventPar_t {
      long long npot;
      int n_good_clusters;
      int n_good_lines;
      int n_good_cosmic_seeds;
      int n_good_time_clusters;
      double weight;
      double gen_energy;

      EventPar_t() {
        init();
      }

      void init(long long np = 0) {
        npot = np;
        n_good_clusters = 0;
        n_good_lines = 0;
        n_good_cosmic_seeds = 0;
        n_good_time_clusters = 0;
        weight = 1.;
        gen_energy = 0.;
      }
    };

    //--------------------------------------------------------------------------------------
    struct SimPar_t {
      const SimParticle* sim;
      unsigned nhits;
      double edep;

      SimPar_t() {
        init();
      }

      void init(const SimParticle* s = nullptr, unsigned n = 0, double e = 0.) {
        sim = s;
        nhits = n;
        edep = e;
      }
    };

    //--------------------------------------------------------------------------------------
    struct ClusterPar_t {
      const CaloCluster*     cluster;
      const KalSeed*         line;
      const CosmicTrackSeed* cosmic_seed;
      const TimeCluster*     time_cluster;
      const CaloClusterMC*   mc;
      const SimParticle*     primary_sim;
      const SimParticle*     secondary_sim;
      const Calorimeter*     calorimeter;

      float r;
      double second_moment;
      double e1;
      double e2;
      double e9;
      double e25;
      double t_var;
      double photon_id;
      int nmatched_lines;
      int nfit_matched_lines;
      int nmatched_cosmic_seeds;
      int nfit_matched_cosmic_seeds;
      int nmatched_time_clusters;
      int nfit_matched_time_clusters;

      float mc_time;
      float sim_1_edep;
      float sim_1_time;
      float sim_1_x;
      float sim_1_y;
      int   sim_1_main_crystal;
      float sim_1_main_crystal_energy;
      int   sim_1_nhits;
      int   sim_1_type;
      int   sim_1_pdg;
      float sim_2_edep;
      float sim_2_time;
      float sim_2_x;
      float sim_2_y;
      int   sim_2_main_crystal;
      float sim_2_main_crystal_energy;
      int   sim_2_nhits;
      int   sim_2_type;
      int   sim_2_pdg;

      ClusterPar_t() {
        init();
        calorimeter = nullptr;
      }

      void init(const CaloCluster* cl = nullptr,  const KalSeed* ln = nullptr,
                const CosmicTrackSeed* cs = nullptr, const TimeCluster* tc = nullptr) {
        cluster = cl;
        line = ln;
        cosmic_seed = cs;
        time_cluster = tc;
        mc = nullptr;
        primary_sim = nullptr;
        secondary_sim = nullptr;
        r = 0.f;
        second_moment = 0.;
        e1 = 0.;
        e2 = 0.;
        e9 = 0.;
        e25 = 0.;
        t_var = 0.;
        photon_id = 0.;
        nmatched_lines = 0;
        nfit_matched_lines = 0;
        nmatched_cosmic_seeds = 0;
        nfit_matched_cosmic_seeds = 0;
        nmatched_time_clusters = 0;
        nfit_matched_time_clusters = 0;
        mc_time = 0.f;
        sim_1_edep = 0.f;
        sim_1_time = 0.f;
        sim_1_x = 0.f;
        sim_1_y = 0.f;
        sim_1_main_crystal = -1;
        sim_1_main_crystal_energy = 0.f;
        sim_1_nhits = 0;
        sim_1_type = -1;
        sim_1_pdg = 0;
        sim_2_edep = 0.f;
        sim_2_time = 0.f;
        sim_2_x = 0.f;
        sim_2_y = 0.f;
        sim_2_main_crystal = -1;
        sim_2_main_crystal_energy = 0.f;
        sim_2_nhits = 0;
        sim_2_type = -1;
        sim_2_pdg = 0;
        if(!cl) return;

        const float x = cluster->cog3Vector().x();
        const float y = cluster->cog3Vector().y();
        r = std::sqrt(x*x + y*y);
        if(calorimeter) {
          ClusterUtils util(*calorimeter, *cluster);
          second_moment = util.secondMoment();
          e1            = util.e1();
          e2            = util.e2();
          e9            = util.e9();
          e25           = util.e25();
        }
        if(cluster->caloHitsPtrVector().size() > 0) {
          const double t = cluster->time();
          for(const auto& hit : cluster->caloHitsPtrVector()) {
            t_var += std::pow(hit->time() - t, 2);
          }
          t_var /= cluster->caloHitsPtrVector().size();
        }
      }

      size_t ncr() const {
        return (cluster) ? cluster->caloHitsPtrVector().size() : size_t(0);
      }

      double frac_1() const {
        if(!cluster) return 0.;
        const double edep = cluster->energyDep();
        if(edep <= 0. || ncr() == 0) return 0.;
        return (cluster->caloHitsPtrVector().at(0)->energyDep() / edep);
      }

      double frac_2() const {
        if(!cluster) return 0.;
        const double edep = cluster->energyDep();
        if(edep <= 0. || ncr() == 0) return 0.;
        if(ncr() == 1) return frac_1();
        return ((cluster->caloHitsPtrVector().at(0)->energyDep() +
                 cluster->caloHitsPtrVector().at(1)->energyDep()) / edep);
      }

      double time_var() const {
        if(!cluster) return -9999.;
        const double t0 = cluster->time();
        const auto& hits = cluster->caloHitsPtrVector();
        if(hits.empty()) return 0.;
        double var = 0.;
        for(const auto& hit : hits) {
          const double dt = hit->time() - t0;
          var += dt*dt;
        }
        var /= hits.size();
        return var;
      }
    };

    //--------------------------------------------------------------------------------------
    struct LinePar_t {
      const KalSeed* line;
      const CosmicTrackSeed* cosmic_seed;
      const TimeCluster* time_cluster;

      LinePar_t() { init(); }
      void init(const KalSeed* l = nullptr) {
        line = l;
        cosmic_seed = nullptr;
        time_cluster = nullptr;
        if(!l) return;
      }
    };

    //--------------------------------------------------------------------------------------
    struct CosmicSeedPar_t {
      const CosmicTrackSeed* seed;

      CosmicSeedPar_t() { init(); }
      void init(const CosmicTrackSeed* s = nullptr) {
        seed = s;
        if(!s) return;
      }
    };

    //--------------------------------------------------------------------------------------
    struct TimeClusterPar_t {
      const TimeCluster* tc;

      int n_primary_hits;
      int n_other_hits;
      int n_total_primary_hits;
      int n_hits_high_z;
      int n_primary_hits_high_z;

      TimeClusterPar_t() { init(); }
      void init(const TimeCluster* t = nullptr) {
        tc = t;
        n_primary_hits = 0;
        n_other_hits = 0;
        n_total_primary_hits = 0;
        n_hits_high_z = 0;
        n_primary_hits_high_z = 0;
        if(!t) return;
      }

      double purity() const {
        if(n_primary_hits + n_other_hits == 0) return 0.;
        return static_cast<double>(n_primary_hits) / (n_primary_hits + n_other_hits);
      }

      double efficiency() const {
        if(n_total_primary_hits == 0) return 0.;
        return static_cast<double>(n_primary_hits) / n_total_primary_hits;
      }
    };

    //--------------------------------------------------------------------------------------
    struct ChargedTrack_t {
      const CaloCluster* cluster;
      const KalSeed* line;

      ChargedTrack_t() { init(); }
      void init(const CaloCluster* c = nullptr, const KalSeed* l = nullptr) {
        cluster = c;
        line = l;
        if(!c || !l) return;
      }
    };

    //--------------------------------------------------------------------------------------
    struct Photon_t {
      const CaloCluster* cluster;
      const KalSeed* line;
      const CosmicTrackSeed* cosmic_seed;
      const TimeCluster* time_cluster;

      float r;

      Photon_t() { init(); }
      void init(const CaloCluster* c = nullptr, const KalSeed* l = nullptr, const CosmicTrackSeed* cs = nullptr, const TimeCluster* tc = nullptr) {
        cluster = c;
        line = l;
        cosmic_seed = cs;
        time_cluster = tc;
        if(!c) return;

        const float x = cluster->cog3Vector().x();
        const float y = cluster->cog3Vector().y();
        r = std::sqrt(x*x + y*y);
      }
    };

} // namespace Run1BAna
