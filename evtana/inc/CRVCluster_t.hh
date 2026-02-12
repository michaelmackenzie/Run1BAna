//
// CRV cluster information
// Michael MacKenzie (2025)

#ifndef RUN1BEVTANA_CRVCLUSTER_T_HH
#define RUN1BEVTANA_CRVCLUSTER_T_HH

// ROOT includes
#include "Rtypes.h"

// EventNtuple includes
#include "EventNtuple/inc/CrvHitInfoReco.hh"
#include "EventNtuple/inc/CrvHitInfoMC.hh"
#include "EventNtuple/inc/RootVectors.hh"

namespace Run1BEvtAna {
  struct CRVCluster_t {
    mu2e::CrvHitInfoReco* crvHit_;
    mu2e::CrvHitInfoMC* crvHitMC_;
    // Note: CRV positions here are translated into track system units
    // Offsets: x_trk = -3904, y_trk = 0, z_trk = 10175

    //-------------------------------------------------
    // Accessors

    int        SectorType()   { return (crvHit_) ? crvHit_->sectorType : 0           ; }
    XYZVectorF Position()     { return (crvHit_) ? crvHit_->pos + XYZVectorF(-3904.f, 0.f, 10175.f) : XYZVectorF(); }
    float      TimeStart()    { return (crvHit_) ? crvHit_->timeStart  : 0.f         ; }
    float      TimeEnd()      { return (crvHit_) ? crvHit_->timeEnd    : 0.f         ; }
    float      Time()         { return (crvHit_) ? crvHit_-> time      : 0.f         ; }
    float      PEs()          { return (crvHit_) ? crvHit_->PEs        : 0.f         ; }
    int        NHits()        { return (crvHit_) ? crvHit_->nHits      : 0           ; }
    int        NLayers()      { return (crvHit_) ? crvHit_->nLayers    : 0           ; }
    float      Slope()        { return (crvHit_) ? crvHit_->angle      : 0.f         ; }

    //-------------------------------------------------
    // Additional functions

    float      PEsPerLayer()  {
      const int nl = NLayers();
      return (nl > 0) ? PEs() / nl : 0.f;
    }
    float      PEsPerHit()  {
      const int nh = NHits();
      return (nh > 0) ? PEs() / nh : 0.f;
    }

    // Extrapolate along the trajectory to a point in z
    float      ZExtrapolation(const float y = 0.f) {
      const float slope = Slope(); // assume dy/dz
      const float y_1(Position().y()), z_1(Position().z());
      if(slope == 0.f) return z_1; // no slope
      // z = (y_1 - y_0) * dz/dy + z_0; y_1 = 0 = solenoid axis
      const float z = (y - y_1) / slope + z_1;
      return z;
    }

    // Distance to given positions
    float Distance(XYZVectorF end) {
      XYZVectorF dist(Position() - end);
      return dist.r();
    }
    float DistanceSTFront() {
      const static float x_st(-3904.f), y_st(0.f), z_st(5471.f); //stopping target entrance position (assume the center of it)
      return Distance(XYZVectorF(x_st, y_st, z_st));
    }
    float DistanceSTBack() {
      const static float x_st(-3904.f), y_st(0.f), z_st(6271.f); //stopping target exit position (assume the center of it)
      return Distance(XYZVectorF(x_st, y_st, z_st));
    }
    float DistanceCaloFront() {
      const static float x_cal(-3904.f), y_cal(0.f), z_cal(11820.f); //calo disk 0 front on solenoid axis
      return Distance(XYZVectorF(x_cal, y_cal, z_cal));
    }
    float DistanceCaloBack() {
      const static float x_cal(-3904.f), y_cal(0.f), z_cal(13220.f); //calo disk 1 back on solenoid axis
      return Distance(XYZVectorF(x_cal, y_cal, z_cal));
    }
    float DistanceExtrapolation(const double y = 0.f, const double x = -3904.f) {
      return Distance(XYZVectorF(x, y, ZExtrapolation(y)));
    }

    // Time to propagate between two points (assuming speed of light)
    float ToF(const float distance) {
      const static float vlightinv = 1.f / 300.f; // light travels ~300 mm / ns
      return vlightinv * distance;
    }
    float TimeAtSTFront  () { return Time() + ToF(DistanceSTFront  ()); }
    float TimeAtSTBack   () { return Time() + ToF(DistanceSTBack   ()); }
    float TimeAtCaloFront() { return Time() + ToF(DistanceCaloFront()); }
    float TimeAtCaloBack () { return Time() + ToF(DistanceCaloBack ()); }
    float TimeAtExtrapolation(const double y = 0.f, const double x = -3904.f) { return Time() + ToF(DistanceExtrapolation(y,x)); }

    // Very approximate times at the tracker front assuming given paths
    float TimeViaSTFront  () { return TimeAtSTFront() + ((Position().z() > 6500.f) ? + 30.f : 10.f); }
    float TimeViaSTBack   () { return TimeAtSTBack () + ((Position().z() > 6500.f) ? + 30.f : 10.f); }
    float TimeViaCaloFront() { return TimeAtCaloFront() + 70.f; }
    float TimeViaCaloBack () { return TimeAtCaloFront() + 80.f; }
    float TimeViaExtrapolation(const double y = 0.f, const double x = -3904.f) {
      const static float vlightinv = 1.f / 300.f; // light travels ~300 mm / ns
      const float z(ZExtrapolation(y));
      const float time_front = std::fabs(z - 8540.f) * vlightinv;
      return TimeAtExtrapolation(y,x) + time_front + ((z > 8500.f) ? 30.f : 10.f);
    }

    void Reset() {
      crvHit_ = nullptr;
    }

    CRVCluster_t() { Reset(); }
  };
}
#endif
