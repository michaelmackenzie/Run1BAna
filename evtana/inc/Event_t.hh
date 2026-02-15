//
// Event information
// Michael MacKenzie (2025)

#ifndef RUN1BEVTANA_EVENT_T_HH
#define RUN1BEVTANA_EVENT_T_HH

// ROOT includes
#include "Rtypes.h"

namespace Run1BEvtAna {
  struct Event_t {
    // Event numbers
    Int_t         run_              ;
    Int_t         subrun_           ;
    Int_t         event_            ;

    // Event info
    Double_t      weight_           ;
    Int_t         ndigis_           ;
    Int_t         ngoodtrks_        ;
    Int_t         ngoodlines_       ;
    Int_t         ntrks_id_         ;
    Int_t         nclusters_        ;
    Float_t       inst_lum_         ;
    Float_t       pb_time_          ;
    Int_t         napr_helices_     ;
    Int_t         ncpr_helices_     ;
    Int_t         noffline_helices_ ;
    Int_t         napr_tracks_      ; // for trigger analysis
    Int_t         ncpr_tracks_      ; // for trigger analysis
    Int_t         ngood_tracks_     ;
    Int_t         ncrv_clusters_    ;
    Int_t         ngood_crvclusters_;

    // Track array counters
    Int_t         ntracks_          ;
    Int_t         nde_tracks_       ;
    Int_t         nue_tracks_       ;
    Int_t         ndmu_tracks_      ;
    Int_t         numu_tracks_      ;
    Int_t         nelectrons_       ;
    Int_t         nmuons_           ;
    Int_t         nprotons_         ;

    Int_t         nlines_           ;

    // Event selection IDs
    Bool_t        passed_apr_       ;
    Bool_t        passed_cpr_       ;
    Bool_t        triggered_        ;
    Int_t         noncrv_vetoid_    ;

    // MC info
    Int_t         nsimps_           ;

    void Reset() {
      weight_            = 1.;
      ndigis_            = 0;
      ngoodtrks_         = 0;
      ntrks_id_          = 0;
      nclusters_         = 0;
      inst_lum_          = 1.f;
      pb_time_           = 0.f;
      napr_helices_      = 0;
      ncpr_helices_      = 0;
      noffline_helices_  = 0;
      napr_tracks_       = 0;
      ncpr_tracks_       = 0;
      ngood_tracks_      = 0;
      ngoodlines_        = 0;
      nde_tracks_        = 0;
      nue_tracks_        = 0;
      numu_tracks_       = 0;
      ndmu_tracks_       = 0;
      ncrv_clusters_     = 0;
      ngood_crvclusters_ = 0;
      ntracks_           = 0;
      nelectrons_        = 0;
      nmuons_            = 0;
      nprotons_          = 0;
      nlines_            = 0;
      passed_apr_        = false;
      passed_cpr_        = false;
      triggered_         = false;
      noncrv_vetoid_     = 0;
      nsimps_            = 0;
    }

    Event_t() { Reset(); }
  };
}
#endif
