#include "TRandom3.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
// generate parameters to misalign the tracker to include an overall twist (dphi/zdz), skew in x and y, and random mis-positions and angles
// David Brown (LBNL)
//
//

using namespace std;
void MisalignTracker(bool misalign_planes, bool misalign_panels, bool misalign_straws, double twist, double skew, double zsqueeze, double rsqueeze, double plsigalpha, double plsigpos, double pasigalpha, double pasigpos, double strawsigv, double strawsigw, const char* outfile, int seed=3499803) {
  unsigned NPlanes(36);
  unsigned NPanels(6);
  unsigned NStraws(96);
  double planezgap(83.0); // approximate number
  double planer(700.0); // approximate number
  double trackerlen = NPlanes*planezgap; // approximate

  TRandom3 myrand(seed);
  double twist_val = myrand.Uniform(-twist,twist);
  double xskew_val = myrand.Uniform(-skew,skew);
  double yskew_val = myrand.Uniform(-skew,skew);
  double zsqueeze_val = myrand.Uniform(-zsqueeze,zsqueeze);
  double rsqueeze_val = myrand.Uniform(-rsqueeze,rsqueeze);
  double dtwist_val = twist_val/(TMath::Pi()*2.0*planer*trackerlen);

  filebuf fb;
  fb.open (outfile,ios::out);
  ostream os(&fb);
  // write header
  os << "# Tracker Misalignments from MisalignTracker.C, parameters: "
  << " seed = "<< seed
  << " twist = "<< twist << " mm" 
  << " skew = "<< skew << " mm" 
  << " zsqueeze = "<< zsqueeze << " mm" 
  << " rsqueeze = "<< rsqueeze << " mm"
  << " planesigalpha = "<< plsigalpha << " rad" 
  << " planesigpos = "<< plsigpos << " mm"
  << " panelsigalpha = "<< pasigalpha << " rad"
  << " panelsigpos = "<< pasigpos << " mm"
  << " strawsigV = "<< strawsigv << " mm"
  << " strawsigW = "<< strawsigw << " mm"
  << " twist_val = "<< twist_val << " mm" 
  << " dtwist_val = " << dtwist_val << " rad/mm"
  << " xskew_val = "<< xskew_val << " mm" 
  << " yskew_val = "<< yskew_val << " mm" 
  << " zsqueeze_val = "<< zsqueeze_val << " mm" 
  << " rsqueeze_val = "<< rsqueeze_val << " mm" 
  << endl;
  os << "# .X MisalignTracker.C(" << misalign_planes << ", " << misalign_panels << ", " << misalign_straws << ", " << twist << ", " << skew << ", " << zsqueeze << ", " << rsqueeze << ", " << plsigalpha << ", " << plsigpos << ", " << pasigalpha << ", " << pasigpos << ", " << strawsigv << ", " << strawsigw <<  ", \"" << outfile << "\", " << seed <<")" << std::endl;
  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;

  // global track: dummy for now
  os << " TABLE TrkAlignTracker " << endl;
  os << "0, 0_0_0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << endl;
  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;
 
  os << "TABLE TrkAlignPlane" << endl;
  for(unsigned iplane=0;iplane < NPlanes; ++iplane){
    if (misalign_planes){
      double planez = (iplane-float(NPlanes)/2.0)*planezgap;
      double dalphax = myrand.Gaus(0.0,plsigalpha);
      double dalphay = myrand.Gaus(0.0,plsigalpha);
      double dalphaz = myrand.Gaus(dtwist_val*planez,plsigalpha);
      double dx = myrand.Gaus(planez*xskew_val/trackerlen,plsigpos);
      double dy = myrand.Gaus(planez*yskew_val/trackerlen,plsigpos);
      double dz = myrand.Gaus(planez*zsqueeze_val/trackerlen,plsigpos);
      if ((iplane % 4) == 1 || (iplane % 4 == 2)){
        dx *= -1;
        dz *= -1;
        dalphax *= -1;
        dalphaz *= -1;
      }

      os << iplane << ", " << iplane << "_0_0, "
        << dx << ", "
        << dy << ", "
        << dz << ", "
        << dalphax << ", "
        << dalphay << ", "
        << dalphaz << endl;
    }else{
      os << iplane << ", " << iplane << "_0_0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << std::endl;
    }
  }
  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;

  os << "TABLE TrkAlignPanel" << endl;
  for(unsigned iplane=0;iplane < NPlanes; ++iplane){
    for(unsigned ipanel=0;ipanel < NPanels; ++ipanel){

      if (misalign_panels){
        double dalphax = myrand.Gaus(0.0,pasigalpha);
        double dalphay = myrand.Gaus(0.0,pasigalpha);
        double dalphaz = myrand.Gaus(0.0,pasigalpha);
        double dx = 0.0;
        double dy = myrand.Gaus(rsqueeze_val,pasigpos);
        double dz = myrand.Gaus(0.0,pasigpos);

        os << iplane*NPanels + ipanel << ", " << iplane << "_" << ipanel << "_0, "
          << dx << ", "
          << dy << ", "
          << dz << ", "
          << dalphax << ", "
          << dalphay << ", "
          << dalphaz << endl;
      }else{
        os << iplane*NPanels + ipanel << ", " << iplane << "_" << ipanel << "_0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0" << std::endl;
      }
    }
  }
  for(unsigned ibuf=0;ibuf<3;ibuf++)
    os << "#" << endl;
  os << "TABLE TrkAlignStraw" << endl;
  for(unsigned iplane=0;iplane < NPlanes; ++iplane){
    for(unsigned ipanel=0;ipanel < NPanels; ++ipanel){
      for(unsigned istraw=0;istraw < NStraws; ++istraw){
        if (misalign_straws){
          double sdvhv = myrand.Gaus(0.0,strawsigv);
          double sdwhv = myrand.Gaus(0.0,strawsigw);
          double wdvhv = myrand.Gaus(0.0,strawsigv);
          double wdwhv = myrand.Gaus(0.0,strawsigw);
          double sdvcal = myrand.Gaus(0.0,strawsigv);
          double sdwcal = myrand.Gaus(0.0,strawsigw);
          double wdvcal = myrand.Gaus(0.0,strawsigv);
          double wdwcal = myrand.Gaus(0.0,strawsigw);
          os << iplane*NPanels*NStraws + ipanel*NStraws + istraw << ", " << iplane << "_" << ipanel << "_" << istraw << ", "
            << wdvcal << ", " << wdwcal << ", "
            << wdvhv << ", " << wdwhv << ", "
            << sdvcal << ", " << sdwcal << ", "
            << sdvhv << ", " << sdwhv << std::endl;
        }else{
          os << iplane*NPanels*NStraws + ipanel*NStraws + istraw << ", " << iplane << "_" << ipanel << "_" << istraw << ", "
            << "0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0" << std::endl;
        }
      }
    }
  }

}
