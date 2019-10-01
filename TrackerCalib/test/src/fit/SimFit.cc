#include <iostream>

#include <TH1F.h>
#include <TRandom.h>

#include "prototype/Utils.hh"

#include "fit/SimFit.hh"

void SimFit::calculate_weighted_pdf (const std::vector<double> &x, TH1F* h, double doca_min, double doca_max, bool dist) const
{
  // x[0] is the tau
  // x[1] is the gaussian width
  // x[2] is time offset between PMT time and straw hit time, so time[i] = x[2] when the ion is directly on the wire
  // x[3] is the background
  double tau = x[0];
  double sigma = x[1];
  double toffset = x[2];
  double background = x[3];

  for (int i=0;i<this->times.size();i++){
    if (this->docas[i] < doca_min || this->docas[i] > doca_max)
      continue;

    double hypotenuse = sqrt(pow(this->docas[i],2) + pow(tau * AVG_VELOCITY,2));
    double tau_eff = hypotenuse/AVG_VELOCITY - this->docas[i]/AVG_VELOCITY;

    for (int j=0;j<h->GetNbinsX();j++){
      double time_residual = h->GetXaxis()->GetBinCenter(j+1);
      if (dist){
        time_residual = h->GetXaxis()->GetBinCenter(j+1)/AVG_VELOCITY;
      }

      double pdf_val = this->interpolatePDF(time_residual,sigma,tau_eff);

      if (pdf_val < 1e-5)
        pdf_val = 1e-5;
      h->SetBinContent(j+1,h->GetBinContent(j+1) + pdf_val);
    }
  }
}



double SimFit::operator() (const std::vector<double> &x) const
{
  // x[0] is the tau
  // x[1] is the gaussian width
  // x[2] is time offset between PMT time and straw hit time, so time[i] = x[2] when the ion is directly on the wire
  // x[3] is the background
  double tau = x[0];
  double sigma = x[1];
  double toffset = x[2];
  double background = x[3];

  if (sigma > pdf_maxs)
    return 1e10;

  double llike = 0;
  for (int i=0;i<this->docas.size();i++){

    double time_residual = this->TimeResidual(this->docas[i],this->times[i],toffset);

    double hypotenuse = sqrt(pow(this->docas[i],2) + pow(tau * AVG_VELOCITY,2));
    double tau_eff = hypotenuse/AVG_VELOCITY - this->docas[i]/AVG_VELOCITY;

//    if (time_residual > 40)
//      time_residual = 40.0;
    if (time_residual < -40)
      time_residual = -40.0;

    //double pdf_val = 1/sqrt(2*TMath::Pi()*x[5]*x[5]) * exp(-(time_residual*time_residual)/(2*x[5]*x[5])) + fabs(x[6]);
    // the above analytic form has lots of floating point issues so we go with a precomputed 3d table and interpolate
    double pdf_val = this->interpolatePDF(time_residual,sigma,tau_eff);

    if (pdf_val < 1e-5)
      pdf_val = 1e-5;

    llike -= log(pdf_val);

  }

  for (int i=0;i<this->nparams;i++){
    if (this->constraints[i] > 0)
      llike += pow((x[i]-this->constraint_means[i])/this->constraints[i],2);
  }

  return llike;
}




