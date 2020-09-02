// Author: S Middleton
// Date: March 2019
// Purpose: Holds functions for the fitting of Cosmic Tracks in tracker

// Mu2e Cosmics:
#include "CosmicReco/inc/CosmicTrackFit.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"

// art
#include "canvas/Persistency/Common/Ptr.h"

//Mu2e General:
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "TrackerGeom/inc/Tracker.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"

//For Drift:
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BbrGeom/Trajectory.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"

//Fitting
#include "Mu2eUtilities/inc/ParametricFit.hh"
#include "Mu2eUtilities/inc/BuildLinearFitMatrixSums.hh"
#include "CosmicReco/inc/MinuitDriftFitter.hh"
#include "CosmicReco/inc/DriftFitUtils.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

//ROOT:
#include "TMatrixD.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "Math/Math.h"
#include "Math/DistFunc.h"

// boost
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

//String:
#include <string>

using namespace std;
using namespace boost::accumulators;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;
using namespace ROOT::Math;

std::vector<double> TrackerConstraints(1500, 1500);

struct ycomp : public std::binary_function<XYZVec,XYZVec,bool> {
    bool operator()(XYZVec const& p1, XYZVec p2) { return p1.y() > p2.y(); }
  };
  

namespace mu2e
{
    CosmicTrackFit::CosmicTrackFit(const Config& conf) :
	_dontuseflag (conf.dontuseflag()),
	_minnsh   (conf.minnsh()),
	_maxniter (conf.maxniter()),
	_maxd (conf.maxd()),
        _useTSeedDirection (conf.UseTSeedDirection())
	{}

    void CosmicTrackFit::Chi2Fit(const char* title, CosmicTrackSeed &tseed, art::Event const& event, ComboHitCollection const& chcol, std::vector<StrawHitIndex> &panelHitIdxs){ 
      // create combohitcollection ordered in y
      ComboHitCollection combohits;
      combohits.setParent(chcol.parent());
      for (size_t i=0;i<panelHitIdxs.size();i++){
        combohits.push_back(chcol[panelHitIdxs[i]]);
      }

      //FIXME
      tseed._status.merge(TrkFitFlag::helixOK);  
      tseed._status.merge(TrkFitFlag::helixConverged);
      //FIXME
      tseed._track.converged = true;



      BuildLinearFitMatrixSums S;
      TrackAxes Axes;
      double a0 = 0, a1 = 0, b0 = 0, b1 = 0;

      // Get Initial Estimate of track direction
      XYZVec currentTrackDirection = (combohits[combohits.size()-1].pos() - combohits[0].pos()).Unit();
      Axes = ParametricFit::GetTrackAxes(currentTrackDirection);

      // Begin iteration for finding the best track fit possible.
      unsigned niter(0);
      bool converged = false;

      // iterate fit until convergence
      while(niter < _maxniter && converged==false){ 
        niter +=1; 

        S.clear(); 
        Axes = ParametricFit::GetTrackAxes(currentTrackDirection);
        for (size_t i=0; i<combohits.size(); ++i){  
          if (!use_hit(combohits[i]))
            continue;
          std::vector<double> ErrorsXY = ParametricFit::GetErrors(combohits[i], Axes._XDoublePrime, Axes._YDoublePrime); 
          S.addPoint(combohits[i].pos(), Axes._XDoublePrime, Axes._YDoublePrime, Axes._ZPrime, ErrorsXY[0], ErrorsXY[1]); 
        }    

        a0 = S.GetAlphaX()[0][0];
        a1 = S.GetAlphaX()[1][0];
        b0 = S.GetAlphaY()[0][0];
        b1 = S.GetAlphaY()[1][0];


        // update cosmic track
        tseed._track.SetFitCoordSystem(Axes);
        TrackParams FitParams(a0,a1,b0,b1);
        tseed._track.SetFitParams(FitParams);

        TrackEquation FitEquation = ConvertFitToDetectorFrame(tseed._track.FitCoordSystem, tseed._track.FitParams.Position(), tseed._track.FitParams.Direction(), true, false);
        tseed._track.SetFitEquation(FitEquation);

        currentTrackDirection = FitEquation.Dir;
        /*
        std::cout << "     alpha: " << S.GetAlphaX()[0][0] << std::endl;
        std::cout << "            " << S.GetAlphaX()[1][0] << std::endl;
        std::cout << "     gamma: " << S.GetGammaX()[0][0] << " " << S.GetGammaX()[0][1] << std::endl;
        std::cout << "            " << S.GetGammaX()[1][0] << " " << S.GetGammaX()[1][1] << std::endl;
        std::cout << "      beta: " << S.GetBetaX()[0][0] << std::endl;
        std::cout << "            " << S.GetBetaX()[1][0] << std::endl;
        std::cout << std::endl;
        std::cout << "     alpha: " << S.GetAlphaY()[0][0] << std::endl;
        std::cout << "            " << S.GetAlphaY()[1][0] << std::endl;
        std::cout << "     gamma: " << S.GetGammaY()[0][0] << " " << S.GetGammaY()[0][1] << std::endl;
        std::cout << "            " << S.GetGammaY()[1][0] << " " << S.GetGammaY()[1][1] << std::endl;
        std::cout << "      beta: " << S.GetBetaY()[0][0] << std::endl;
        std::cout << "            " << S.GetBetaY()[1][0] << std::endl;
        */

        // Identify outliers
        size_t worsti = 0;
        double worstchi2 = 0;
        for (size_t i=0; i<combohits.size(); ++i){  
          if (!use_hit(combohits[i]))
            continue;
          std::vector<double> ErrorsXY = ParametricFit::GetErrors(combohits[i], Axes._XDoublePrime, Axes._YDoublePrime); 
          double hitchi2 = ParametricFit::GetHitChi2(a0, a1, ErrorsXY[0], combohits[i].pos(), Axes._XDoublePrime, Axes._ZPrime) + ParametricFit::GetHitChi2(b0, b1, ErrorsXY[1], combohits[i].pos(), Axes._YDoublePrime, Axes._ZPrime);
          if (hitchi2 > worstchi2){
            worstchi2 = hitchi2;
            worsti = i;
          }
        }
        if (worstchi2 > 9){ //FIXME
          combohits[worsti]._flag.merge(StrawHitFlag::outlier);
        }else{
          break;
        }
      }


      // make sure _straw_chits is the straw level hits not panel hits
      // FIXME any reason to do here?
      //std::vector<ComboHitCollection::const_iterator> chids;  
      //chcol.fillComboHits(event, panelHitIdxs, chids); 
      //for (auto const& it : chids){
      //  tseed._straw_chits.push_back(it[0]);
      //}
      tseed._straw_chits = combohits;
      return;
    }


	 
   

  
    XYZVec CosmicTrackFit::ConvertPointToDetFrame(XYZVec vec){
        Hep3Vector vec1(vec.x(),vec.y(),vec.z());
        GeomHandle<DetectorSystem> det;
        Hep3Vector vec2 = det->toDetector(vec1);
	XYZVec XYZ(vec2.x(), vec2.y(), vec2.z());
	return XYZ;

    }

/*------------Translate fit back into XYZ and/or the detector frame------//
Using matrices to ctransform from local to global coordinates
//-----------------------------------------------------------------------*/
TrackEquation CosmicTrackFit::ConvertFitToDetectorFrame(TrackAxes axes, XYZVec Position, XYZVec Direction, bool seed, bool det){
	TMatrixD A(3,3);
	A[0][0] = axes._XDoublePrime.X();
	A[0][1] = axes._YDoublePrime.X();
	A[0][2] = axes._ZPrime.X();

	A[1][0] = axes._XDoublePrime.Y();
	A[1][1] = axes._YDoublePrime.Y();
	A[1][2] = axes._ZPrime.Y();

	A[2][0] = axes._XDoublePrime.Z();
	A[2][1] = axes._YDoublePrime.Z();
	A[2][2] = axes._ZPrime.Z();

	TMatrixD P(3,1);
	P[0][0] = Position.X();
	P[1][0] = Position.Y();
	P[2][0] = Position.Z();

	TMatrixD D(3,1);
	D[0][0] = Direction.X();
	D[1][0] = Direction.Y();
	D[2][0] = Direction.Z();

	TMatrixD PXYZ(A*P);
	TMatrixD DXYZ(A*D);
	
        XYZVec Pos(PXYZ[0][0], PXYZ[1][0], PXYZ[2][0]);
	XYZVec Dir(DXYZ[0][0], DXYZ[1][0] , DXYZ[2][0]);

        Pos -= Dir * Pos.y()/Dir.y();

    //    std::cout << "Converting: " << axes._XDoublePrime << " " << axes._YDoublePrime << " " << axes._ZPrime << std::endl;
    //    std::cout << "  " << Position << " " << Direction << std::endl;
    //    std::cout << "  " << Pos << " " << Dir << std::endl;
    //    std::cout << std::endl;

        if (det == true){ // is this detector frame?
          XYZVec PosInDet = ConvertPointToDetFrame(Pos);
          return TrackEquation(PosInDet, Dir);
        } else {
          return TrackEquation(Pos, Dir);
        }
}


    bool CosmicTrackFit::goodTrack(CosmicTrack& track) 
    { 
      return true;
//	if(track.Diag.FinalChiTot < _max_seed_chi2) return true;
//	else return false; 
    }

    bool CosmicTrackFit::use_hit(const ComboHit& thit) const 
    {
        return (!thit._flag.hasAnyProperty(_dontuseflag) && thit.nStrawHits() >= _minnsh);
    }

    bool CosmicTrackFit::use_track(double track_length) const 
    {
	return (track_length > _maxd) ? false : true ;
    }

}//end namespace
