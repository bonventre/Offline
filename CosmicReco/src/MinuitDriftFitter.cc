// Author : S Middleton
// Date : August 2019
// Purpose: calls  minuit fitting to cosmic track seed. Input is CosmicTrackSeed, can then derive
// parameters from CosmicTrack stored there.

// ROOT:
#include "CosmicReco/inc/MinuitDriftFitter.hh"
#include "CosmicReco/inc/PDFFit.hh"
#include "Math/Math.h"
#include "Math/VectorUtil.h"
#include "Mu2eUtilities/inc/ParametricFit.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "TMath.h"
#include "TrackerGeom/inc/Tracker.hh"

// For Drift:
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/BbrGeom/Trajectory.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
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

void RunMigrad(
    std::vector<double> & pars, 
    std::vector<double> & errors,
    std::vector<double> & cov_out,
    bool & minuit_converged,
    ROOT::Minuit2::FCNBase const& fit,
    int diag, double mntolerance, double mnprecision) {
  
  // Initiate Minuit Fit:
  ROOT::Minuit2::MnStrategy mnStrategy(2);
  ROOT::Minuit2::MnUserParameters params(pars, errors);
  ROOT::Minuit2::MnMigrad migrad(fit, params, mnStrategy);

  if (mnprecision > 0) {
    migrad.SetPrecision(mnprecision);
  }

  // Define Minimization method as "MIGRAD" (see minuit documentation)
  ROOT::Minuit2::FunctionMinimum min = migrad(0, mntolerance);
  if (diag > 1) {
    ROOT::Minuit2::MnPrint::SetLevel(3);
    ROOT::Minuit2::operator<<(cout, min);
  } else {
    ROOT::Minuit2::MnPrint::SetLevel(0);
  }

  // Will be the results of the fit routine:
  ROOT::Minuit2::MnUserParameters const& results = min.UserParameters();

  minuit_converged = min.IsValid();
  pars = results.Params();
  errors = results.Errors();

  cov_out = min.UserCovariance().Data();
}

void DoDriftTimeFit(int const& diag, CosmicTrackSeed& tseed, StrawResponse const& srep,
                    const Tracker* tracker, double mntolerance, double mnprecision) {

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
  RunMigrad(pars, errors, tseed._track.MinuitParams.cov, 
    tseed._track.minuit_converged, fit, 
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

  XYZVec X(1, 0, 0);
  XYZVec Y(0, 1, 0);
  XYZVec Z(0, 0, 1);

  TrackAxes XYZ(X, Y, Z);
  tseed._track.MinuitCoordSystem = XYZ;
  tseed._track.MinuitEquation.Pos =
      XYZVec(tseed._track.MinuitParams.A0, 0, tseed._track.MinuitParams.B0);
  tseed._track.MinuitEquation.Dir =
      XYZVec(tseed._track.MinuitParams.A1, -1, tseed._track.MinuitParams.B1);

  for (size_t i = 0; i < tseed._straw_chits.size(); i++) {
    Straw const& straw = tracker->getStraw(tseed._straw_chits[i].strawId());
    TwoLinePCA pca(straw.getMidPoint(), straw.getDirection(),
                   Geom::Hep3Vec(tseed._track.MinuitEquation.Pos),
                   Geom::Hep3Vec(tseed._track.MinuitEquation.Dir));
    if (pca.dca() > 2.5) {
      tseed._straw_chits[i]._flag.merge(StrawHitFlag::outlier);
    }
  }
}

void DoSimpleDriftFit(int const& diag, CosmicTrackSeed& tseed, StrawResponse const& srep,
                    const Tracker* tracker, double simpleTres, double mntolerance, double mnprecision) {

  auto dir = tseed._track.FitEquation.Dir;
  auto intercept = tseed._track.FitEquation.Pos;
  dir /= -1 * dir.y();
  intercept -= dir * intercept.y() / dir.y();

  // now gaussian fit, transverse distance only
  std::vector<double> errors(4, 0);
  std::vector<double> pars(4, 0);

  pars[0] = intercept.x();
  pars[1] = intercept.z();
  pars[2] = dir.x();
  pars[3] = dir.z();
  errors[0] = tseed._track.FitParams.Covarience.sigA0;
  errors[1] = tseed._track.FitParams.Covarience.sigB0;
  errors[2] = tseed._track.FitParams.Covarience.sigA1;
  errors[3] = tseed._track.FitParams.Covarience.sigB1;

  // Define the PDF used by Minuit:
  SimpleDriftFit fit(tseed._straw_chits, srep, tracker, simpleTres);
  RunMigrad(pars, errors, tseed._track.MinuitParams.cov, 
    tseed._track.minuit_converged, fit, 
    diag, mntolerance, mnprecision);

  tseed._track.MinuitParams.A0 = pars[0];
  tseed._track.MinuitParams.B0 = pars[1];
  tseed._track.MinuitParams.A1 = pars[2];
  tseed._track.MinuitParams.B1 = pars[3];
  tseed._track.MinuitParams.T0 = fit.t0(pars);
  tseed._track.MinuitParams.deltaA0 = errors[0];
  tseed._track.MinuitParams.deltaB0 = errors[1];
  tseed._track.MinuitParams.deltaA1 = errors[2];
  tseed._track.MinuitParams.deltaB1 = errors[3];
  tseed._track.MinuitParams.deltaT0 = tseed._t0._t0err;
  tseed._t0._t0 = tseed._track.MinuitParams.T0;
  tseed._t0._t0err = tseed._track.MinuitParams.deltaT0;

  XYZVec X(1, 0, 0);
  XYZVec Y(0, 1, 0);
  XYZVec Z(0, 0, 1);

  TrackAxes XYZ(X, Y, Z);
  tseed._track.MinuitCoordSystem = XYZ;
  tseed._track.MinuitEquation.Pos =
      XYZVec(tseed._track.MinuitParams.A0, 0, tseed._track.MinuitParams.B0);
  tseed._track.MinuitEquation.Dir =
      XYZVec(tseed._track.MinuitParams.A1, -1, tseed._track.MinuitParams.B1);

  for (size_t i = 0; i < tseed._straw_chits.size(); i++) {
    Straw const& straw = tracker->getStraw(tseed._straw_chits[i].strawId());
    TwoLinePCA pca(straw.getMidPoint(), straw.getDirection(),
                   Geom::Hep3Vec(tseed._track.MinuitEquation.Pos),
                   Geom::Hep3Vec(tseed._track.MinuitEquation.Dir));
    if (pca.dca() > 2.5) {
      tseed._straw_chits[i]._flag.merge(StrawHitFlag::outlier);
    }
  }
}

} // namespace MinuitDriftFitter
