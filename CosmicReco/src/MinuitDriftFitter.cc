// Author : S Middleton
// Date : August 2019
// Purpose: calls  minuit fitting to cosmic track seed. Input is CosmicTrackSeed, can then derive
// parameters from CosmicTrack stored there.

// ROOT:
#include "Offline/CosmicReco/inc/MinuitDriftFitter.hh"
#include "Offline/CosmicReco/inc/PDFFit.hh"
#include "Math/Math.h"
#include "Math/VectorUtil.h"
#include "Offline/Mu2eUtilities/inc/ParametricFit.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "TMath.h"
#include "Offline/TrackerGeom/inc/Tracker.hh"

#include <TObjString.h>
#include <TROOT.h>
#include <TSystem.h>

// Minuit
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>

using namespace mu2e;

namespace MinuitDriftFitter {

void DoDriftTimeFit(
    std::vector<double> & pars,
    std::vector<double> & errors,
    std::vector<double> & cov_out,
    bool & minuit_converged,
    GaussianDriftFit &fit,
    double driftres,
    int diag, double mntolerance, double mnprecision) {

  // Initiate Minuit Fit:
  ROOT::Minuit2::MnStrategy mnStrategy(2);
  ROOT::Minuit2::MnUserParameters params(pars, errors);
  ROOT::Minuit2::MnMigrad migrad(fit, params, mnStrategy);

  if (mnprecision > 0) {
    migrad.SetPrecision(mnprecision);
  }

  // Do first fit stage with fixed drift res
  // and minimal t0
  fit.setFixedT0(true);
  fit.setFixedDriftRes(true,driftres);
  migrad.Fix(4);
  ROOT::Minuit2::FunctionMinimum temp_min = migrad(0, mntolerance);
  ROOT::Minuit2::MnUserParameters const& temp_results = temp_min.UserParameters();
  for (size_t i=0;i<4;i++){
    pars[i] = temp_results.Params()[i];
    errors[i] = temp_results.Errors()[i];
  }
  XYZVectorF ft0pos(pars[0], 0, pars[1]);
  XYZVectorF ft0dir(pars[2], -1, pars[3]);
  ft0dir = ft0dir.unit();
  pars[4] = fit.averageT0(pars);
  fit.setFixedT0(false);
  fit.setFixedDriftRes(false);
  migrad.Release(4);

  // Define Minimization method as "MIGRAD" (see minuit documentation)
  ROOT::Minuit2::FunctionMinimum min = migrad(0, mntolerance);
  if (diag > 1) {
    ROOT::Minuit2::MnPrint::SetGlobalLevel(3);
    ROOT::Minuit2::operator<<(std::cout, min);
  } else {
    ROOT::Minuit2::MnPrint::SetGlobalLevel(0);
  }

  // Will be the results of the fit routine:
  ROOT::Minuit2::MnUserParameters const& results = min.UserParameters();

  minuit_converged = min.IsValid();
  pars = results.Params();
  errors = results.Errors();

  if (min.HasValidCovariance()) {
    cov_out = min.UserCovariance().Data();
  } else {
    cov_out = std::vector<double>(15, 0);
  }
}

void DoDriftTimeFit(int const& diag, CosmicTrackSeed& tseed, StrawResponse const& srep,
                    const Tracker* tracker, double driftres, double mntolerance, double mnprecision) {

  auto dir = tseed._track.FitEquation.Dir;
  auto intercept = tseed._track.FitEquation.Pos;
  dir /= -1 * dir.y();
  intercept -= dir * intercept.y() / dir.y();

  // now gaussian fit, transverse distance only
  std::vector<double> errors(5, 0);
  std::vector<double> pars(5, 0);

  pars[0] = intercept.x();
  pars[1] = intercept.z();
  pars[2] = dir.x();
  pars[3] = dir.z();
  pars[4] = tseed._t0._t0;
  errors[0] = tseed._track.FitParams.Covarience.sigA0;
  errors[1] = tseed._track.FitParams.Covarience.sigB0;
  errors[2] = tseed._track.FitParams.Covarience.sigA1;
  errors[3] = tseed._track.FitParams.Covarience.sigB1;
  errors[4] = tseed._t0.t0Err();

  // Define the PDF used by Minuit:
  GaussianDriftFit fit(tseed._straw_chits, srep, tracker);
  DoDriftTimeFit(pars, errors, tseed._track.MinuitParams.cov,
    tseed._track.minuit_converged, fit, driftres,
    diag, mntolerance, mnprecision);

  tseed._track.MinuitParams.A0 = pars[0];
  tseed._track.MinuitParams.B0 = pars[1];
  tseed._track.MinuitParams.A1 = pars[2];
  tseed._track.MinuitParams.B1 = pars[3];
  tseed._track.MinuitParams.T0 = pars[4];
  tseed._track.MinuitParams.deltaA0 = errors[0];
  tseed._track.MinuitParams.deltaB0 = errors[1];
  tseed._track.MinuitParams.deltaA1 = errors[2];
  tseed._track.MinuitParams.deltaB1 = errors[3];
  tseed._track.MinuitParams.deltaT0 = errors[4];
  tseed._t0._t0 = tseed._track.MinuitParams.T0;
  tseed._t0._t0err = tseed._track.MinuitParams.deltaT0;

  XYZVectorF X(1, 0, 0);
  XYZVectorF Y(0, 1, 0);
  XYZVectorF Z(0, 0, 1);

  TrackAxes XYZ(X, Y, Z);
  tseed._track.MinuitCoordSystem = XYZ;
  tseed._track.MinuitEquation.Pos =
      XYZVectorF(tseed._track.MinuitParams.A0, 0, tseed._track.MinuitParams.B0);
  tseed._track.MinuitEquation.Dir =
      XYZVectorF(tseed._track.MinuitParams.A1, -1, tseed._track.MinuitParams.B1);

  for (size_t i = 0; i < tseed._straw_chits.size(); i++) {
    Straw const& straw = tracker->getStraw(tseed._straw_chits[i].strawId());
    TwoLinePCA pca(straw.getMidPoint(), straw.getDirection(),
                   GenVector::Hep3Vec(tseed._track.MinuitEquation.Pos),
                   GenVector::Hep3Vec(tseed._track.MinuitEquation.Dir));
    if (pca.dca() > 2.5) {
      tseed._straw_chits[i]._flag.merge(StrawHitFlag::outlier);
    }
  }
}

} // namespace MinuitDriftFitter
