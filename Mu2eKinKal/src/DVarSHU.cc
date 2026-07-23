#include "Offline/Mu2eKinKal/inc/DVarSHU.hh"
#include <cmath>

namespace mu2e {
  using KinKal::ClosestApproachData;
  using KinKal::VEC3;
  DVarSHU::DVarSHU(Config const& config) {
    maxdoca_ = std::get<0>(config);
    double maxdocaerr = std::get<1>(config);
    maxdvar_ = maxdocaerr*maxdocaerr;
    minrdrift_ = std::get<2>(config);
    maxrdrift_ = std::get<3>(config);
    minsigma_ = std::get<4>(config);
    maxchi2_ = std::get<5>(config);
    maxchi2full_ = std::get<6>(config);
    std::string flag = std::get<7>(config);
    flag_ = KKSHFlag(flag);
    std::string allowed = std::get<8>(config);
    allowed_ = WHSMask(allowed);
    std::string freeze = std::get<9>(config);
    freeze_ = WHSMask(freeze);
    diag_ = std::get<10>(config);
    if(diag_ > 0)std::cout << "DVarSHU max doca, doca error " << maxdoca_ << " " << maxdocaerr
      << " rdrift range [" << minrdrift_ << "," << maxrdrift_ << "] Allowing "
        << allowed_ << " Freezing " << freeze_ << " Flags " << flag
        << std::endl;
  }

  WireHitState DVarSHU::wireHitState(WireHitState const& input, ClosestApproachData const& tpdata,DriftInfo const& dinfo, double T, double d_meas, double var_d, double d_meas2, double var_d2, double wprior2, double &wplus, double &wminus, double &wplus2, double &wminus2) const {
    WireHitState whstate = input;
    bool updated(false);
    if(input.updateable(StrawHitUpdaters::CAD)){
      whstate.quality_[WireHitState::sign] = 0;

      double absdoca = fabs(tpdata.doca());
      if(dinfo.rDrift_ < maxrdrift_ && dinfo.rDrift_ > minrdrift_ && absdoca < maxdoca_ && tpdata.docaVar() > 0.0 && tpdata.docaVar() < maxdvar_ ){
        T = std::fabs(minsigma_);

        double rplus = d_meas - tpdata.doca();
        double rminus = -1*d_meas - tpdata.doca();
//        if (d_meas < 0){
//          double sgn = (tpdata.doca() >= 0.0) ? 1.0 : -1.0;
//          rplus = -tpdata.doca() + sgn * d_meas;
//          rminus = -tpdata.doca() + sgn * d_meas;
//        }

        double totvar = var_d + tpdata.docaVar();
        double chi2plus = rplus*rplus/totvar;
        double chi2minus = rminus*rminus/totvar;
        double aplus = -1*chi2plus/(2*T);
        double aminus = -1*chi2minus/(2*T);

        if (wprior2 > 0){
          double rplus2 = d_meas2 - tpdata.doca();
          double rminus2 = -1*d_meas2 - tpdata.doca();
//          if (d_meas2 < 0){
//            double sgn = (tpdata.doca() >= 0.0) ? 1.0 : -1.0;
//            rplus2 = -tpdata.doca() + sgn * d_meas2;
//            rminus2 = -tpdata.doca() + sgn * d_meas2;
//          }
          double totvar2 = var_d2 + tpdata.docaVar();
          double chi2plus2 = rplus2*rplus2/totvar;
          double chi2minus2 = rminus2*rminus2/totvar;
          double aplus2 = -1*chi2plus2/(2*T);
          double aminus2 = -1*chi2minus2/(2*T);

          double logwplus = std::log(1-wprior2) - 0.5 * log(totvar) + aplus;
          double logwminus = std::log(1-wprior2) - 0.5 * log(totvar) + aminus;
          double logwplus2 = std::log(wprior2) - 0.5 * log(totvar2) + aplus2;
          double logwminus2 = std::log(wprior2) - 0.5 * log(totvar2) + aminus2;
          double logwmax = std::max({logwplus, logwminus, logwplus2, logwminus2});

          double eplus = std::exp(logwplus - logwmax);
          double eminus = std::exp(logwminus - logwmax);
          double eplus2 = std::exp(logwplus2 - logwmax);
          double eminus2 = std::exp(logwminus2 - logwmax);
          double esum = eplus + eminus + eplus2 + eminus2;

          wplus = eplus / esum;
          wplus2 = eplus2 / esum;
          wminus = eminus / esum;
          wminus2 = eminus2 / esum;
        }else{
          double amax = std::max(aplus,aminus);
          double eplus = std::exp(aplus - amax);
          double eminus = std::exp(aminus - amax);
          wplus = eplus / (eplus + eminus);
          wminus = eminus / (eplus + eminus);
          wplus2 = 0;
          wminus2 = 0;
        }

        whstate.wplus_ = wplus;
        whstate.wplus2_ = wplus2;
        whstate.wminus_ = wminus;
        whstate.wminus2_ = wminus2;

        if ((rplus*rplus/var_d > maxchi2_ && rminus*rminus/var_d > maxchi2_) || (chi2plus > maxchi2full_ && chi2minus > maxchi2full_)){
          whstate.state_ = WireHitState::null;
          updated = true;
        }else{
          if(absdoca/sqrt(tpdata.docaVar()) > minsigma_){
            // in the sweet spot: use the DOCA to sign the ambiguity
            if(allowed_.hasAnyProperty(WHSMask::drift)) {
              whstate.state_ = tpdata.doca() > 0.0 ? WireHitState::right : WireHitState::left;
              updated = true;
            }
          } else if(allowed_.hasAnyProperty(WHSMask::null)) {
            whstate.state_ = WireHitState::null;
            updated = true;
          }
        }
      } else if(allowed_.hasAnyProperty(WHSMask::inactive)) {
        whstate.state_ = WireHitState::inactive;
        updated = true;
      }
      if(updated){
        whstate.algo_ = StrawHitUpdaters::CAD;
        whstate.flag_ = flag_;
        whstate.frozen_ = whstate.isIn(freeze_);
      }
      if (diag_ > 1)std::cout << "DVarSHU set hit " << whstate << std::endl;
    } else if (diag_ > 1) {
      std::cout << "DVarSHU skipping hit " << whstate << std::endl;    }
    return whstate;
  }

  std::string const& DVarSHU::configDescription() {
    static std::string descrip( "Maximum DOCA to use hit, Maximum DOCA error to use hit, Minimum rdrift to set LR ambiguity, Maximum rdrift to use hit, allowed states, States to freeze, diag level");
    return descrip;
  }

}
