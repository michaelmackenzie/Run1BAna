//
// A single kinematic line
// Michael MacKenzie (2026)

#ifndef RUN1BEVTANA_LINE_T_HH
#define RUN1BEVTANA_LINE_T_HH

// ROOT includes
#include "Rtypes.h"
#include "TString.h"

// Event Ntuple includes
#include "EventNtuple/inc/TrkInfo.hh"
#include "EventNtuple/inc/TrkCaloHitInfo.hh"
#include "EventNtuple/inc/TrkSegInfo.hh"
#include "EventNtuple/inc/LoopHelixInfo.hh"
#include "EventNtuple/inc/SimInfo.hh"

#include "EventNtuple/rooutil/inc/Track.hh"

// local includes
#include "Run1BAna/evtana/inc/GlobalConstants.h"
#include "Run1BAna/evtana/inc/CRVCluster_t.hh"
#include "Run1BAna/evtana/inc/CaloCluster_t.hh"

namespace Run1BEvtAna {
  struct Line_t {

    rooutil::Track* line_      = nullptr; // pointer to the EventNtuple::Track object
    CRVCluster_t* stub_        = nullptr; // best matched CRV cluster
    CaloCluster_t* cluster_    = nullptr; // best matched calo cluster

    // Line IDs
    int id_[kMaxTrackIDs];

    // Fit observables
    double obs_[kMaxObservables];

    //----------------------------------------------
    // Line information accessors

    //----------------------------------------------
    // Basic line fit checks
    bool IsGood() {
      if(!line_) return false;
      if(!line_->trk) return false;
      if(line_->trk->status < 0) return false;
      if(line_->trk->goodfit == 0) return false;
      return true;
    }

    //----------------------------------------------
    // Line fit info
    float Chi2Dof() { return (line_ && line_->trk && line_->trk->ndof > 0) ? line_->trk->chisq / line_->trk->ndof : -1.f; }
    float FitCon () { return (line_ && line_->trk) ? line_->trk->fitcon : -1.f; }
    int   NActive() { return (line_ && line_->trk) ? line_->trk->nactive : -1; }
    float TrkQual() { return (line_ && line_->trkqual && line_->trkqual->valid) ? line_->trkqual->result : -1000.f; }

    //----------------------------------------------
    // Particle hypothesis used in the fit
    int FitPDG() {
      if(!line_ || !line_->trk) return 0;
      return line_->trk->pdg;
    }
    int Charge() {
      const int pdg = FitPDG();
      const int abs_pdg = std::abs(pdg);
      const bool negative = std::signbit(pdg); // true if negative
      switch(abs_pdg) {
      case 11: case 13: case 15:
        return (negative) ?  1 : -1;
      case 211: case 2212:
        return (negative) ? -1 :  1;
      case 0: return 0;
      default:
        break;
      }
      printf("Line_t::%s: Unknown fit hypothesis %i --> Returning 0 charge\n", __func__, pdg);
      return 0;
    }

    //----------------------------------------------
    // MC Particle info
    const mu2e::SimInfo* SimInfo() {
      if(!line_) return nullptr;
      if(!line_->trkmcsim) return nullptr;

      // Find the sim particle with the most active hits
      int max_hits(-1);
      const mu2e::SimInfo* sim_info(nullptr);
      for(const auto& info : *(line_->trkmcsim)) {
        if(info.nactive > max_hits) {
          max_hits = info.nactive;
          sim_info = &info;
        }
      }
      return sim_info;
    }

    int MCPDG    () { auto sim_info = SimInfo(); return (sim_info) ? sim_info->pdg            : 0; }
    int MCHits   () { auto sim_info = SimInfo(); return (sim_info) ? sim_info->nhits          : 0; }
    int MCActive () { auto sim_info = SimInfo(); return (sim_info) ? sim_info->nactive        : 0; }
    int MCProcess() { auto sim_info = SimInfo(); return (sim_info) ? sim_info->startCode      : 0; }
    int MCGenP   () { auto sim_info = SimInfo(); return (sim_info) ? sim_info->mom.r()        : 0; }
    int MCGenE   () {
      auto sim_info = SimInfo();
      if(!sim_info) return 0.;
      // Retrieve the particle mass and add it to the momentum (if available)
      const double mass((std::abs(sim_info->pdg) > 10000) ? 0. : ParticleMass(sim_info->pdg));
      return std::sqrt(std::pow(sim_info->mom.r(), 2) + mass*mass);
    }

    //----------------------------------------------
    // Line-Calo hit retrieval
    const mu2e::TrkCaloHitInfo* TCH() { return (line_) ? line_->trkcalohit : nullptr; }

    //----------------------------------------------
    // Segment info
    const mu2e::TrkSegInfo* Segment(mu2e::SurfaceIdDetail::enum_type surface) {
      if(!line_ || !line_->trksegs) return nullptr;
      for(const auto& seg : *(line_->trksegs)) {
        if(seg.sid == surface) return &seg;
      }
      return nullptr;
    }
    int SegmentIndex(mu2e::SurfaceIdDetail::enum_type surface) {
      if(!line_ || !line_->trksegs) return -1;
      for(size_t index = 0; index < line_->trksegs->size(); ++index) {
        const auto& seg =  line_->trksegs->at(index);
        if(seg.sid == surface) return int(index);
      }
      return -1;
    }
    const mu2e::LoopHelixInfo* LHSegment(mu2e::SurfaceIdDetail::enum_type surface) {
      if(!line_ || !line_->trksegpars_lh) return nullptr;
      const int index = SegmentIndex(surface);
      if(index < 0 || index >= int(line_->trksegpars_lh->size())) return nullptr;
      return &(line_->trksegpars_lh->at(index));
    }
    const mu2e::SurfaceStepInfo* MCSegment(mu2e::SurfaceIdDetail::enum_type surface) {
      if(!line_ || !line_->trksegsmc) return nullptr; // No MC info
      auto reco_seg = Segment(surface);
      if(!reco_seg) return nullptr; // Can't do anything if there's no reco info

      // Search for an MC segment matched to this reco segment
      const mu2e::SurfaceStepInfo* seg(nullptr);
      for(const auto& mc_seg : *(line_->trksegsmc)) {
        if(mc_seg.sid != surface) continue;
        // Found this MC surface
        if(seg) { // check which is segment is a better match
          const float dt_curr = std::fabs(reco_seg->time - seg->time);
          const float dt_new  = std::fabs(reco_seg->time - mc_seg.time);
          if(dt_curr > dt_new) seg = &mc_seg;
        } else seg = &mc_seg;
      }
      return seg;
    }

    //----------------------------------------------
    // Lineer front segment info
    const mu2e::TrkSegInfo* FrontSeg() { return Segment(mu2e::SurfaceIdDetail::TT_Front); }

    //----------------------------------------------
    // Line kinematics at a given surface
    float PSegment     (mu2e::SurfaceIdDetail::enum_type surface) { auto seg =   Segment(surface); return (seg) ? seg->mom.r()               :  0.; }
    float PTSegment    (mu2e::SurfaceIdDetail::enum_type surface) { auto seg =   Segment(surface); return (seg) ? seg->mom.rho()             : -1.; }
    float PZSegment    (mu2e::SurfaceIdDetail::enum_type surface) { auto seg =   Segment(surface); return (seg) ? seg->mom.z()               :  0.; }
    float DMomSegment  (mu2e::SurfaceIdDetail::enum_type surface) { auto seg =   Segment(surface); return (seg) ? seg->dmom                  :  0.; }
    float MomErrSegment(mu2e::SurfaceIdDetail::enum_type surface) { auto seg =   Segment(surface); return (seg) ? seg->momerr                : -1.; }
    float TSegment     (mu2e::SurfaceIdDetail::enum_type surface) { auto seg =   Segment(surface); return (seg) ? seg->time                  :  0.; }
    float TErrSegment  (mu2e::SurfaceIdDetail::enum_type surface) { auto seg = LHSegment(surface); return (seg) ? seg->t0err                 : -1.; }
    float RMaxSegment  (mu2e::SurfaceIdDetail::enum_type surface) { auto seg = LHSegment(surface); return (seg) ? seg->maxr                  :  0.; }
    float RadiusSegment(mu2e::SurfaceIdDetail::enum_type surface) { auto seg = LHSegment(surface); return (seg) ? std::fabs(seg->rad)        :  0.; }
    float RMinSegment  (mu2e::SurfaceIdDetail::enum_type surface) { return std::fabs(RMaxSegment(surface) - 2.f*RadiusSegment(surface)); }

    float D0Segment(mu2e::SurfaceIdDetail::enum_type surface) {
      auto seg = LHSegment(surface);
      // return (seg) ? seg->d0 : -1.e6;
      if(!seg) return -1.e6;
      // FIXME: Evaluating this locally to get the sign
      const double radius = std::fabs(seg->rad);
      const double max_r  = seg->maxr;
      return max_r - 2.*radius;
    }

    float TanDipSegment(mu2e::SurfaceIdDetail::enum_type surface) {
      auto seg = LHSegment(surface);
      if(seg) return seg->tanDip;
      // Estimate it by hand if the segment is not available
      // FIXME: Check that this definition is correct
      // tan(theta) = x/y = pz / pt
      const float pt = PTSegment(surface);
      const float pz = PZSegment(surface);
      if(pt <= 0. || pz == 0.) return -100.;
      return pz/pt;
    }

    //----------------------------------------------
    // MC Line kinematics at a given surface
    float MCPSegment     (mu2e::SurfaceIdDetail::enum_type surface) { auto seg = MCSegment(surface); return (seg) ? seg->mom.r()               :  0.; }
    float MCPTSegment    (mu2e::SurfaceIdDetail::enum_type surface) { auto seg = MCSegment(surface); return (seg) ? seg->mom.rho()             : -1.; }
    float MCPZSegment    (mu2e::SurfaceIdDetail::enum_type surface) { auto seg = MCSegment(surface); return (seg) ? seg->mom.z()               :  0.; }
    float MCTSegment     (mu2e::SurfaceIdDetail::enum_type surface) { auto seg = MCSegment(surface); return (seg) ? seg->time                  :  0.; }

    //----------------------------------------------
    // Line kinematics at the lineer front
    float PFront     () { return PSegment     (mu2e::SurfaceIdDetail::TT_Front); }
    float PTFront    () { return PTSegment    (mu2e::SurfaceIdDetail::TT_Front); }
    float PZFront    () { return PZSegment    (mu2e::SurfaceIdDetail::TT_Front); }
    float DMomFront  () { return DMomSegment  (mu2e::SurfaceIdDetail::TT_Front); }
    float MomErrFront() { return MomErrSegment(mu2e::SurfaceIdDetail::TT_Front); }
    float TFront     () { return TSegment     (mu2e::SurfaceIdDetail::TT_Front); }
    float TErrFront  () { return TErrSegment  (mu2e::SurfaceIdDetail::TT_Front); }
    float D0Front    () { return D0Segment    (mu2e::SurfaceIdDetail::TT_Front); }
    float TanDipFront() { return TanDipSegment(mu2e::SurfaceIdDetail::TT_Front); }
    float RMaxFront  () { return RMaxSegment  (mu2e::SurfaceIdDetail::TT_Front); }
    float RadiusFront() { return RadiusSegment(mu2e::SurfaceIdDetail::TT_Front); }
    int   Trajectory () {
      const float pz = PZFront();
      if(pz == 0.f) return 0;
      return (pz < 0.f) ? -1 : 1;
    }

    float MCPFront     () { return MCPSegment   (mu2e::SurfaceIdDetail::TT_Front); }
    float MCPTFront    () { return MCPTSegment  (mu2e::SurfaceIdDetail::TT_Front); }
    float MCPZFront    () { return MCPZSegment  (mu2e::SurfaceIdDetail::TT_Front); }
    float MCTFront     () { return MCTSegment   (mu2e::SurfaceIdDetail::TT_Front); }
    float MCDeltaPFront() { return PFront() - MCPFront(); }
    int   MCTrajectory () {
      const float pz = MCPZFront();
      if(pz == 0.f) return 0;
      return (pz < 0.f) ? -1 : 1;
    }

    //----------------------------------------------
    // Line kinematics at the lineer middle
    float PMiddle     () { return PSegment     (mu2e::SurfaceIdDetail::TT_Mid); }
    float PTMiddle    () { return PTSegment    (mu2e::SurfaceIdDetail::TT_Mid); }
    float PZMiddle    () { return PZSegment    (mu2e::SurfaceIdDetail::TT_Mid); }
    float DMomMiddle  () { return DMomSegment  (mu2e::SurfaceIdDetail::TT_Mid); }
    float MomErrMiddle() { return MomErrSegment(mu2e::SurfaceIdDetail::TT_Mid); }
    float TMiddle     () { return TSegment     (mu2e::SurfaceIdDetail::TT_Mid); }
    float TErrMiddle  () { return TErrSegment  (mu2e::SurfaceIdDetail::TT_Mid); }
    float D0Middle    () { return D0Segment    (mu2e::SurfaceIdDetail::TT_Mid); }
    float TanDipMiddle() { return TanDipSegment(mu2e::SurfaceIdDetail::TT_Mid); }
    float RMaxMiddle  () { return RMaxSegment  (mu2e::SurfaceIdDetail::TT_Mid); }

    //----------------------------------------------
    // Line kinematics at the lineer back
    float PBack     () { return PSegment     (mu2e::SurfaceIdDetail::TT_Back); }
    float PTBack    () { return PTSegment    (mu2e::SurfaceIdDetail::TT_Back); }
    float PZBack    () { return PZSegment    (mu2e::SurfaceIdDetail::TT_Back); }
    float DMomBack  () { return DMomSegment  (mu2e::SurfaceIdDetail::TT_Back); }
    float MomErrBack() { return MomErrSegment(mu2e::SurfaceIdDetail::TT_Back); }
    float TBack     () { return TSegment     (mu2e::SurfaceIdDetail::TT_Back); }
    float TErrBack  () { return TErrSegment  (mu2e::SurfaceIdDetail::TT_Back); }
    float D0Back    () { return D0Segment    (mu2e::SurfaceIdDetail::TT_Back); }
    float TanDipBack() { return TanDipSegment(mu2e::SurfaceIdDetail::TT_Back); }
    float RMaxBack  () { return RMaxSegment  (mu2e::SurfaceIdDetail::TT_Back); }

    //----------------------------------------------
    // Line kinematics at the stopping target exit
    float PSTBack     () { return PSegment     (mu2e::SurfaceIdDetail::ST_Back); }
    float PTSTBack    () { return PTSegment    (mu2e::SurfaceIdDetail::ST_Back); }
    float PZSTBack    () { return PZSegment    (mu2e::SurfaceIdDetail::ST_Back); }
    float DMomSTBack  () { return DMomSegment  (mu2e::SurfaceIdDetail::ST_Back); }
    float MomErrSTBack() { return MomErrSegment(mu2e::SurfaceIdDetail::ST_Back); }
    float TSTBack     () { return TSegment     (mu2e::SurfaceIdDetail::ST_Back); }
    float TErrSTBack  () { return TErrSegment  (mu2e::SurfaceIdDetail::ST_Back); }
    float D0STBack    () { return D0Segment    (mu2e::SurfaceIdDetail::ST_Back); }
    float TanDipSTBack() { return TanDipSegment(mu2e::SurfaceIdDetail::ST_Back); }
    float RMaxSTBack  () { return RMaxSegment  (mu2e::SurfaceIdDetail::ST_Back); }

    float MCPSTBack     () { return MCPSegment   (mu2e::SurfaceIdDetail::ST_Back); }
    float MCPTSTBack    () { return MCPTSegment  (mu2e::SurfaceIdDetail::ST_Back); }
    float MCPZSTBack    () { return MCPZSegment  (mu2e::SurfaceIdDetail::ST_Back); }
    float MCTSTBack     () { return MCTSegment   (mu2e::SurfaceIdDetail::ST_Back); }
    float MCDeltaPSTBack() { return PSTBack() - MCPSTBack(); }

    //----------------------------------------------
    // Line CaloHit info
    float ECluster() { auto tch = TCH(); return (!tch) ?    0. : tch->edep  ; }
    float Dt      () { auto tch = TCH(); return (!tch) ? -1.e6 : tch->dt    ; }
    float DOCA    () { auto tch = TCH(); return (!tch) ? -1.e6 : tch->doca  ; }
    float PTOCA   () { auto tch = TCH(); return (!tch) ? -1.e6 : tch->ptoca ; }
    float CDepth  () { auto tch = TCH(); return (!tch) ? -1.e6 : tch->cdepth; }
    float EPFront () {
      const float ecl(ECluster()), p(PFront());
      return (p > 0.) ? ecl/p : 0.;
    }

    //----------------------------------------------
    // Accessing/setting the line IDs
    void SetID(const int ID, const int index = 0) {
      if(index < 0 || index >= kMaxTrackIDs) throw std::runtime_error(Form("Accessing a line ID index (%i) out of bounds!", index));
      id_[index] = ID;
    }
    int ID(const int index = 0) {
      if(index < 0 || index >= kMaxTrackIDs) throw std::runtime_error(Form("Accessing line ID index (%i) out of bounds!", index));
      return id_[index];
    }

    //----------------------------------------------
    // Accessing/setting the observable
    void SetObs(const double val, const int index = 0) {
      if(index < 0 || index >= kMaxObservables) throw std::runtime_error(Form("Accessing an observable index (%i) out of bounds!", index));
      obs_[index] = val;
    }
    double Obs(const int index = 0) {
      if(index < 0 || index >= kMaxObservables) throw std::runtime_error(Form("Accessing an observable index (%i) out of bounds!", index));
      return obs_[index];
    }

    //----------------------------------------------
    // Reset the input info
    void Reset() {
      line_    = nullptr;
      stub_    = nullptr;
      cluster_ = nullptr;
      for(int iid = 0; iid < kMaxTrackIDs; ++iid) id_[iid] = 0;
      for(int iobs = 0; iobs < kMaxObservables; ++iobs) obs_[iobs] = 0.;
    }

    //----------------------------------------------
    // Print the line
    void Print(TString opt = "") {
      opt.ToLower();
      if(opt.Contains("banner")) {
        std::string filler(130, '-');
        printf("%s\n", filler.c_str());
        printf("Idx: %5s %10s %10s %10s %10s %7s %6s %10s %5s %5s %8s %8s\n", "Hyp", "p", "pT", "pz", "t", "Ecl", "tandip", "p(MC)", "PDG", "good", "fitcon", "trkqual");
        printf("%s\n", filler.c_str());
      }
      if(!line_) return;
      printf("Idx: %5i %10.2f %10.2f %10.2f %10.1f %7.1f %6.2f %10.1f %5i %5i %.2e %8.5f\n", FitPDG(), PFront(), PTFront(), PZFront(), TFront(), ECluster(), TanDipFront(),
             MCPFront(), MCPDG(), IsGood(), FitCon(), TrkQual());
    }

    Line_t() { Reset(); }
  };
}
#endif
