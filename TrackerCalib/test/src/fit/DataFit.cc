#include <iostream>
#include <numeric>

#include <TH1F.h>
#include <TRandom.h>
#include <TMath.h>

#include "prototype/Utils.hh"

#include "fit/DataFit.hh"

double DataFit::operator() (const std::vector<double> &x) const
{
  // x[0] is the tau
  // x[1] is the gaussian width
  // x[2] is time offset between PMT time and straw hit time, so time[i] = x[2] when the ion is directly on the wire
  // x[3] is the background
  // x[4] is distance from back of pixel to wire horizontally perpendicular to the straw in mm
  // x[5] is distance from bot pixel to wire in mm
  // x[6] is the pixel separation
  double tau = x[0];
  double sigma = x[1];
  double toffset = x[2];
  double background = x[3];
  double hdist = x[4];
  double vdist = x[5];
  double pixel_separation = x[6];

  if (sigma > pdf_maxs)
    return 1e10;

  double llike = 0;
  for (int i=0;i<this->times.size();i++){
    double doca = calculate_DOCA(this->tops[i],this->bots[i],hdist,vdist,pixel_separation);
    double time_residual = this->TimeResidual(doca,this->times[i],toffset);
    double doca_penalty = 1;
    
    if (doca >= 2.5){
      doca_penalty = exp(-(doca-2.5)/5.);
      doca = 2.5;
    }


    double hypotenuse = sqrt(pow(doca,2) + pow(tau * AVG_VELOCITY,2));
    double tau_eff = hypotenuse/AVG_VELOCITY - doca/AVG_VELOCITY;

//    if (time_residual > 40)
//      time_residual = 40.0;
    if (time_residual < -40)
      time_residual = -40.0;

    //double pdf_val = 1/sqrt(2*TMath::Pi()*x[5]*x[5]) * exp(-(time_residual*time_residual)/(2*x[5]*x[5])) + fabs(x[6]);
    // the above analytic form has lots of floating point issues so we go with a precomputed 3d table and interpolate
    double pdf_val = this->interpolatePDF(time_residual,sigma,tau_eff);

    pdf_val *= doca_penalty;

    if (pdf_val < 1e-5)
      pdf_val = 1e-5;

    llike -= log(pdf_val);

  }

  for (int i=0;i<7;i++){
    if (this->constraints[i] > 0)
      llike += pow((x[i]-this->constraint_means[i])/this->constraints[i],2);
  }
//  std::cout << llike << " : " << tau << " " << sigma << " " << toffset << " " << background << " " << hdist << " " << vdist << " " << pixel_separation << " " << times.size() << std::endl;

  return llike;
}


void DataFit::calculate_weighted_pdf (const std::vector<double> &x, TH1F* h, double doca_min, double doca_max, bool dist) const
{
  // x[0] is the tau
  // x[1] is the gaussian width
  // x[2] is time offset between PMT time and straw hit time, so time[i] = x[2] when the ion is directly on the wire
  // x[3] is the background
  // x[4] is distance from back of pixel to wire horizontally perpendicular to the straw in mm
  // x[5] is distance from bot pixel to wire in mm
  // x[6] is the pixel separation
  double tau = x[0];
  double sigma = x[1];
  double toffset = x[2];
  double background = x[3];
  double hdist = x[4];
  double vdist = x[5];
  double pixel_separation = x[6];

  for (int i=0;i<this->times.size();i++){
    double doca = calculate_DOCA(this->tops[i],this->bots[i],hdist,vdist,pixel_separation);
    double doca_penalty = 1;
    
    if (doca >= 2.5){
      doca_penalty = exp(-(doca-2.5)/5.);
      doca = 2.5;
    }

    if (doca < doca_min || doca > doca_max)
      continue;


    double hypotenuse = sqrt(pow(doca,2) + pow(tau * AVG_VELOCITY,2));
    double tau_eff = hypotenuse/AVG_VELOCITY - doca/AVG_VELOCITY;

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

void DataFit::metropolis(std::vector<double> &seed, std::vector<double> &errors, std::vector<std::vector<double> > &results, std::vector<double> &nlls, int steps)
{
  int accepted = 0;
  std::vector<double> current_x = seed;
  std::vector<double> proposed_x = seed;
  std::vector<double> jumps = errors;
  std::vector<std::vector<double> > jump_vectors;
  for (int j=0;j<7;j++){
    jump_vectors.push_back(std::vector<double>());
  }
  double current_nll = this->operator()(current_x);
  for (int i=0;i<steps*0.1;i++){
    for (int j=0;j<7;j++){
      if (j == 3)
        continue;
      proposed_x[j] = current_x[j] + gRandom->Gaus(0,jumps[j]);
    }
    double proposed_nll = this->operator()(proposed_x);
    double u = gRandom->Uniform();
    if (proposed_nll < current_nll || exp(current_nll - proposed_nll) >= u){
      accepted++;
      current_nll = proposed_nll;
      for (int j=0;j<7;j++){
        current_x[j] = proposed_x[j];
      }
    }
    for (int j=0;j<7;j++){
      jump_vectors[j].push_back(current_x[j]);
    }
  }
  std::cout << "During burn in accepted " << accepted/(steps*0.1) << std::endl;
  for (int j=0;j<7;j++){
    double sum = std::accumulate(jump_vectors[j].begin(), jump_vectors[j].end(), 0.0);
    double mean = sum / jump_vectors[j].size();
    double sq_sum = std::inner_product(jump_vectors[j].begin(), jump_vectors[j].end(), jump_vectors[j].begin(), 0.0);
    double stdev = std::sqrt(sq_sum / jump_vectors[j].size() - mean * mean);
    jumps[j] = stdev *2.4*2.4/7.;
    std::cout << "Retuned jump " << j << " to " << jumps[j] << std::endl;
    jump_vectors[j] = std::vector<double>();
  }
  accepted = 0;
  for (int i=0;i<steps*0.1;i++){
    for (int j=0;j<7;j++){
      if (j == 3)
        continue;
      proposed_x[j] = current_x[j] + gRandom->Gaus(0,jumps[j]);
    }
    double proposed_nll = this->operator()(proposed_x);
    double u = gRandom->Uniform();
    if (proposed_nll < current_nll || exp(current_nll - proposed_nll) >= u){
      accepted++;
      current_nll = proposed_nll;
      for (int j=0;j<7;j++){
        current_x[j] = proposed_x[j];
      }
    }
    for (int j=0;j<7;j++){
      jump_vectors[j].push_back(current_x[j]);
    }
  }
  std::cout << "During burn in accepted " << accepted/(steps*0.1) << std::endl;
  for (int j=0;j<7;j++){
    double sum = std::accumulate(jump_vectors[j].begin(), jump_vectors[j].end(), 0.0);
    double mean = sum / jump_vectors[j].size();
    double sq_sum = std::inner_product(jump_vectors[j].begin(), jump_vectors[j].end(), jump_vectors[j].begin(), 0.0);
    double stdev = std::sqrt(sq_sum / jump_vectors[j].size() - mean * mean);
    jumps[j] = stdev *2.4*2.4/7./2.;
    std::cout << "Retuned jump " << j << " to " << jumps[j] << std::endl;
    jump_vectors[j] = std::vector<double>();
  }

  double minnll = 1e10;
  std::vector<double> bestfit;
  
  accepted = 0;
  for (int i=0;i<steps;i++){
    if (i%10000 == 0)
      std::cout << i << " / " << steps << std::endl;
    for (int j=0;j<7;j++){
      if (j == 3)
        continue;
      proposed_x[j] = current_x[j] + gRandom->Gaus(0,jumps[j]);
    }
    double proposed_nll = this->operator()(proposed_x);
    double u = gRandom->Uniform();
    if (proposed_nll < current_nll || exp(current_nll - proposed_nll) >= u){
      accepted++;
      current_nll = proposed_nll;
      for (int j=0;j<7;j++){
        current_x[j] = proposed_x[j];
      }
    }
    if (current_nll < minnll){
      minnll = current_nll;
      bestfit = current_x;
  //    std::cout << "minnll: " << minnll << " : ";
  //    for (int j=0;j<7;j++){
  //      std::cout << bestfit[j] << " ";
  //    }
  //    std::cout << std::endl;
    }
    results.push_back(current_x);
    nlls.push_back(current_nll);
  }
  std::cout << "Final acceptance: " << accepted/(steps+1e-8) << std::endl;
  this->operator()(bestfit);
}

double GaussianDataFit::operator() (const std::vector<double> &x) const
{
  // x[0] is the tau (NOT USED)
  // x[1] is the gaussian width
  // x[2] is time offset between PMT time and straw hit time, so time[i] = x[2] when the ion is directly on the wire
  // x[3] is the background
  // x[4] is distance from back of pixel to wire horizontally perpendicular to the straw in mm
  // x[5] is distance from bot pixel to wire in mm
  // x[6] is the pixel separation
  double tau = x[0];
  double sigma = x[1];
  double toffset = x[2];
  double background = x[3];
  double hdist = x[4];
  double vdist = x[5];
  double pixel_separation = x[6];

  double llike = 0;
  for (int i=0;i<this->tops.size();i++){
    if (this->outlier[i])
      continue;

    double doca = calculate_DOCA(this->tops[i],this->bots[i],hdist,vdist,pixel_separation);
    double time_residual = this->TimeResidual(doca,this->times[i],toffset);
    
    if (time_residual > 40)
      time_residual = 40.0;
    if (time_residual < -40)
      time_residual = -40.0;

    double pdf_val = 1/sqrt(2*TMath::Pi()*sigma*sigma) * exp(-(time_residual*time_residual)/(2*sigma*sigma));

    llike -= log(pdf_val);
    
  }

  for (int i=0;i<7;i++){
    if (this->constraints[i] > 0)
      llike += pow((x[i]-this->constraint_means[i])/this->constraints[i],2);
  }
  return llike;
}

void GaussianDataFit::SelectOutliers(std::vector<double> &x)
{
  double tau = x[0];
  double sigma = x[1];
  double toffset = x[2];
  double background = x[3];
  double hdist = x[4];
  double vdist = x[5];
  double pixel_separation = x[6];

  for (int i=0;i<this->tops.size();i++){
    double doca = calculate_DOCA(this->tops[i],this->bots[i],hdist,vdist,pixel_separation);
    double time_residual = this->TimeResidual(doca,this->times[i],toffset);
    if (fabs(time_residual) > 10)
      this->outlier[i] = true;
    else
      this->outlier[i] = false;
  }
}
 
