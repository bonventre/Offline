#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>

#include <TCanvas.h>
#include <TMath.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TF1.h>
#include <TFile.h>
#include <TLine.h>
#include <TString.h>

#include "prototype/Utils.hh"
#include "fit/StrawFit.hh"

double StrawFit::operator() (const std::vector<double> &x) const
{
  // x[0] is mean efficiency in straw
  // x[1] is exponential outside straw
  // x[2] is distance from back of pixel to center of straw horizontally perpendicular to the straw in mm
  // x[3] is distance from bot pixel to center of straw in mm

  double meaneff = x[0];
  double exponent = x[1];
  double hdist = x[2];
  double vdist = x[3];

  double llike = 0;
  for (int i=0;i<this->alltops.size();i++){
    if (this->has_shower[i])
      continue;

    double dist = calculate_DOCA(this->alltops[i],this->allbots[i],hdist,vdist);
//    if (dist > 4 && this->has_straw[i]){
//      continue;
//    }
    double pdf_prob = fabs(meaneff) * (1 - 0.5*(1+TMath::Erf((dist-STRAW_RADIUS)/(fabs(exponent)*sqrt(2)))));
    if (pdf_prob > 1.0)
      pdf_prob = 0.999;

    if (pdf_prob < 1e-3 || dist > 3.0)
      pdf_prob = 1e-3;

    if (this->has_straw[i]){
      llike -= log(pdf_prob); 
    }else{
      llike -= log(1-pdf_prob);
    }
  }
//  if (x[2] < 0)
//    llike += 100*x[2]*x[2];
//  if (x[3] < 0)
//    llike += 100*x[3]*x[3];

//  std::cout << llike << " :  " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl; 
 return llike;
}

double StrawFitPL::operator() (const std::vector<double> &x) const
{
  // x[0] is distance from back of pixel to center of straw horizontally perpendicular to the straw in mm
  // x[1] is distance from bot pixel to center of straw in mm
  // x[2] is mean efficiency in straw
  // x[3] is path length where equal to threshold on average
  // x[4] is gaussian width to threshold


  double llike = 0;
  for (int i=0;i<this->alltops.size();i++){
    if (this->has_shower[i])
      continue;

    double dist = calculate_DOCA(this->alltops[i],this->allbots[i],x[0],x[1]);
    double path_length = sqrt(2.485*2.485-dist*dist)*2;
    if (dist > 2.485)
      path_length = -1*(dist-2.485);
//    if (dist > 4 && this->has_straw[i]){
//      continue;
//    }
    
    double pdf_prob = x[2] * (0.5*(1+TMath::Erf((path_length-x[3])/(fabs(x[4])*sqrt(2))))) + 1e-3;
    //std::cout << i << " " << dist << " " << path_length << " " << (path_length-x[3])/(fabs(x[4])*sqrt(2)) << " " << pdf_prob << std::endl;
    if (this->has_straw[i]){
      llike -= log(pdf_prob); 
    }else{
      llike -= log(1-pdf_prob);
    }
  }
  if (x[2] < 0)
    llike += 100*x[2]*x[2];
  if (x[3] < 0)
    llike += 100*x[3]*x[3];

//  std::cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " : " << llike << std::endl; 
  return llike;
}

double StrawFitPL2::operator() (const std::vector<double> &x) const
{
  // x[0] is distance from back of pixel to center of straw horizontally perpendicular to the straw in mm
  // x[1] is distance from bot pixel to center of straw in mm
  // x[2] is mean efficiency in straw
  // x[3] is polya stretching


  double llike = 0;
  for (int i=0;i<this->alltops.size();i++){
    if (this->has_shower[i])
      continue;

    double dist = calculate_DOCA(this->alltops[i],this->allbots[i],x[0],x[1]);
    double path_length = sqrt(2.485*2.485-dist*dist)*2;
    if (dist > 2.485){
      double pdf_prob = 1e-2; //*(1/(1+(dist-2.485)));
      if (this->has_straw[i]){
        llike -= log(pdf_prob);
      }else{
        llike -= log(1-pdf_prob);
      }
      continue;
    }

    double pdf_prob;
    if (path_length*x[3] >= 10)
      pdf_prob = x[2];
    else
      pdf_prob = x[2] * hpcdf->Interpolate(path_length*x[3]) + 1e-2;
    
    if (this->has_straw[i]){
      llike -= log(pdf_prob); 
    }else{
      llike -= log(1-pdf_prob);
    }
  }

//  std::cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << " : " << llike << std::endl; 
  return llike;
}
