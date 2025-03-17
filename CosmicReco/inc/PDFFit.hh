#ifndef _COSMIC_RECO_PDFFit_HH
#define _COSMIC_RECO_PDFFit_HH

#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"

// Tracker Details:
#include "Offline/TrackerConditions/inc/StrawDrift.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

// ROOT
#include "TF1.h"
#include "TH1F.h"
#include "TMath.h"

// Minuit
#include <Minuit2/FCNBase.h>

using namespace mu2e;

class GaussianDriftFit : public ROOT::Minuit2::FCNBase {
public:
  ComboHitCollection shs;
  StrawResponse const& srep;
  const Tracker* tracker;

  int excludeHit;

  bool fixedT0;
  bool fixedDriftRes;
  double driftRes;
  bool constrainToStraw;

  GaussianDriftFit(ComboHitCollection const& _shs, StrawResponse const& _srep,
                   const Tracker* _tracker) :
      shs(_shs),
      srep(_srep), tracker(_tracker), excludeHit(-1),
      fixedT0(false), fixedDriftRes(false), driftRes(0), constrainToStraw(true) {};
  // this tells Minuit to scale variances as if operator() returns a chi2 instead of a log
  // likelihood
  double Up() const { return 1.0; };
  double operator()(const std::vector<double>& x) const;

  void setExcludeHit(int const& hitIdx) {
    excludeHit = hitIdx;
  }
  void setFixedT0(bool fix) { fixedT0 = fix;}
  void setFixedDriftRes(bool fix, double dr=10) { fixedDriftRes = fix; driftRes = dr; }
  void setConstrainToStraw(bool constrain) { constrainToStraw = constrain; }

  double averageT0(const std::vector<double> &x) const;

  double DOCAresidual(ComboHit const& sh, CosmicTrackSeed const& tseed) const {
    std::vector<double> x = {tseed._track.MinuitParams.A0, tseed._track.MinuitParams.B0,
                             tseed._track.MinuitParams.A1, tseed._track.MinuitParams.B1,
                             tseed._track.MinuitParams.T0};
    return DOCAresidual(sh, x);
  }
  double TimeResidual(ComboHit const& sh, CosmicTrackSeed const& tseed) const {
    std::vector<double> x = {tseed._track.MinuitParams.A0, tseed._track.MinuitParams.B0,
                             tseed._track.MinuitParams.A1, tseed._track.MinuitParams.B1,
                             tseed._track.MinuitParams.T0};
    return TimeResidual(sh, x);
  }
  double DOCAresidualError(ComboHit const& sh, CosmicTrackSeed const& tseed) const {
    std::vector<double> x = {tseed._track.MinuitParams.A0, tseed._track.MinuitParams.B0,
                             tseed._track.MinuitParams.A1, tseed._track.MinuitParams.B1,
                             tseed._track.MinuitParams.T0};
    return DOCAresidualError(sh, x, tseed._track.MinuitParams.cov);
  }

  int HitAmbiguity(ComboHit const& sh, CosmicTrackSeed const& tseed) const {
    std::vector<double> x {tseed._track.MinuitParams.A0, tseed._track.MinuitParams.B0,
                             tseed._track.MinuitParams.A1, tseed._track.MinuitParams.B1,
                             tseed._track.MinuitParams.T0};
    return HitAmbiguity(sh, x);
  }

  double reduced_chisq(const std::vector<double>& x);

  int HitAmbiguity(ComboHit const& sh, const std::vector<double>& x) const;

  double DOCAresidual(ComboHit const& sh, const std::vector<double>& x) const;

  double TimeResidual(ComboHit const& sh, const std::vector<double>& x) const;

  double DOCAresidualError(ComboHit const& sh, const std::vector<double>& x,
                           const std::vector<double>& cov) const;
};

#endif
