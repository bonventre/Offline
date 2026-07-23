#ifndef Mu2eKinKal_DVarSHU_hh
#define Mu2eKinKal_DVarSHU_hh
//
// Simple updater of StrawHits based on box cuts of Closest Approach (CA) and drift information
//
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include "Offline/TrackerConditions/inc/DriftInfo.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include "Offline/Mu2eKinKal/inc/KKSHFlag.hh"
#include <tuple>
#include <string>
#include <iostream>

namespace mu2e {
  // Update based just on PTCA to the wire
  class DVarSHU {
    public:
      using Config = std::tuple<float,float,float,float,float,float,float,std::string,std::string,std::string,int>;
      DVarSHU(Config const& config);
      static std::string const& configDescription(); // description of the variables
      // set the state based on the current PTCA value
      WireHitState wireHitState(WireHitState const& input, KinKal::ClosestApproachData const& tpdata,DriftInfo const& dinfo, double T, double d_meas, double var_d, double d_meas2, double var_d2, double wprior2, double& wplus, double& wminus, double& wplus2, double& wminus2) const;
    private:
      double maxdoca_ =0; // maximum DOCA to use hit
      double maxdvar_ =0; // maximum DOCA variance to use hit
      double minrdrift_ =0; // minimum rdrift to use drift information
      double maxrdrift_ =0; // maximum rdrift to use hit
      double minsigma_ = 0;
      double maxchi2_ = 0;
      double maxchi2full_ = 0;
      WHSMask allowed_; // allowed states
      WHSMask freeze_; // states to freeze
      KKSHFlag flag_; // flags
      int diag_ =0; // diag print level
  };
}
#endif
