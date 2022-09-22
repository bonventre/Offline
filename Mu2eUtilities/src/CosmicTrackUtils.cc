//
// Helper class to hold functions formerly in RecoDataProducts/inc/CosmicTrack.hh but which
// use classes from Mu2eUtilities.
//
#include "TMath.h"
#include <iostream>

#include "Offline/Mu2eUtilities/inc/CosmicTrackUtils.hh"

#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA_XYZ.hh"

namespace mu2e {

  std::tuple <double, double, double, double, double, double> KinKalTrackParams( CosmicTrack const& trk){
    XYZVectorF zpos(0.,0.,0);
    XYZVectorF zdir(0.,0.,1.);
    XYZVectorF pos0(trk.MinuitParams.A0, 0, trk.MinuitParams.B0);
    XYZVectorF dir(trk.MinuitParams.A1, -1, trk.MinuitParams.B1);
    dir = dir.Unit();

    std::tuple <double,double, double, double, double, double> info;
    TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(pos0, dir, zpos, zdir);
    XYZVectorF POCA = PCA.point1()-PCA.point2();
    double DOCA = PCA.dca();
    double amsign = copysign(1.0, -(zdir.Cross(POCA)).Dot(dir));

    double d0 = amsign*DOCA;
//    double phi0 = dir.Phi();
    double phi0 = PCA.point1().Phi() - TMath::Pi()/2.0;
    double z0 = PCA.point1().Z();
    double cost = dir.Z();
    double t0 = trk.MinuitParams.T0; //TODO
    double mom = 2000.0;//TODO
    info = std::make_tuple(d0,phi0,z0,cost, t0, mom);

    std::cout << "PCA POINT1 " << PCA.point1().x() << " " << PCA.point1().y() << " " << PCA.point1().z() << " : " << PCA.point1().Phi() << std::endl;
    std::cout << "Was setting to " << dir.Phi() << " " << dir.Z() << " , Now setting to " << PCA.point1().Phi() - 3.1415926535/2. << " " << dir.Z() << std::endl;
//    std::cout << pars[KTRAJ::d0_] << " " << pars[KTRAJ::phi0_] << " " << pars[KTRAJ::z0_] << " " << pars[KTRAJ::cost_] << " " << pars[KTRAJ::t0_] << " " << pars[KTRAJ::mom_] << std::en

    return info;
  }

}
