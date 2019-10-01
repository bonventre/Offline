#include <TMath.h>
#include <cmath>

#include "prototype/Utils.hh"

#include "fit/TrackerCalibFit.hh"

TableFit::TableFit(int _k, std::vector<double> &_constraint_means, std::vector<double> &_constraints, double voltage) : Fit(_constraint_means, _constraints, voltage), k(_k)
{
  pdf_times = new double[pdf_tbins];
  pdf_taus = new double[pdf_taubins];
  pdf_sigmas = new double[pdf_sbins];
  pdf = new double[pdf_tbins*pdf_taubins*pdf_sbins];
  pdf_mint = -10;
  pdf_mintau = 0.1;
  pdf_mins = 0.1;
  pdf_maxt = 40;
  pdf_maxtau = 20.;
  pdf_maxs = 10.;
  pdf_deltat = (pdf_maxt-pdf_mint)/((double) pdf_tbins);
  pdf_deltas = (pdf_maxs-pdf_mins)/((double) pdf_sbins);
  pdf_deltatau = (pdf_maxtau-pdf_mintau)/((double) pdf_taubins);
  for (int i=0;i<pdf_tbins;i++)
    pdf_times[i] = pdf_mint + pdf_deltat * i;
  for (int i=0;i<pdf_sbins;i++)
    pdf_sigmas[i] = pdf_mins + pdf_deltas * i;
  for (int i=0;i<pdf_taubins;i++)
    pdf_taus[i] = pdf_mintau + pdf_deltatau * i;

  this->calculate_full_pdf();
};



double Fit::TimeResidual(double doca, double time, double time_offset) const
{
    double expected_time = this->D2T(doca,time_offset);
    return time-expected_time;
}

double Fit::D2T(double distance, double time_offset) const
{
    double expected_time = time_offset;
    expected_time += sd->D2T(distance);
    return expected_time;
}

double Fit::T2D(double time, double time_offset) const
{
  for (int i=0;i<3000;i++){
    double dist = i*0.001;
    double time_pred = this->D2T(dist,time_offset);
    if (time_pred > time)
      return dist;
  }
  return 2.5;
}


int factorial(int k)
{
  if (k == 0)
    return 1;
  int response = 1;
  for (int i=1;i<k;i++)
    response *= i;
  return response;
}

void TableFit::calculate_full_pdf() {
  for (int is=0;is<pdf_sbins;is++){
    double sigma = this->pdf_sigmas[is];
    for (int it0=0;it0<pdf_tbins;it0++){
      double time_gaus = this->pdf_times[it0];
      double val_gaus = 1.0/sqrt(2*TMath::Pi()*sigma*sigma)*exp(-(time_gaus*time_gaus)/(2*sigma*sigma));

      for (int itau=0;itau<pdf_taubins;itau++){
        double tau = this->pdf_taus[itau];
        for (int it1=0;it1<pdf_tbins-it0;it1++){
          double time_tau = this->pdf_deltat*it1;
          double val_tau = pow(1/tau,k)*pow(time_tau,k-1)*exp(-time_tau/tau)/(double) factorial(k-1);
          this->pdf[is * pdf_taubins * pdf_tbins + itau * pdf_tbins + (it0+it1)] += val_gaus * val_tau;
        }
      }
    }
  }

  for (int is=0;is<pdf_sbins;is++){
    for (int itau=0;itau<pdf_taubins;itau++){
      double total = 0;
      for (int it=0;it<pdf_tbins;it++){
        total += this->pdf[is * pdf_taubins * pdf_tbins + itau * pdf_tbins + it];
      }
      for (int it=0;it<pdf_tbins;it++){
        this->pdf[is * pdf_taubins * pdf_tbins + itau * pdf_tbins + it] /= total;
      }
    }
  }
}


double TableFit::interpolatePDF(double time_residual, double sigma, double tau) const
{
  int bin_s = (sigma - this->pdf_mins)/(this->pdf_deltas);
  if (bin_s >= pdf_sbins-1)
    bin_s = pdf_sbins-2;
  if (bin_s < 0)
    bin_s = 0;
  double s_d = (sigma - (bin_s*this->pdf_deltas+this->pdf_mins))/(this->pdf_deltas);

  int bin_tau = (tau - this->pdf_mintau)/(this->pdf_deltatau);
  if (bin_tau >= pdf_taubins-1)
    bin_tau = pdf_taubins-2;
  double tau_d = (tau - (bin_tau*this->pdf_deltatau+this->pdf_mintau))/(this->pdf_deltatau);

  int bin_t = (time_residual - this->pdf_mint)/(this->pdf_deltat);
  if (bin_t >= pdf_tbins-1)
    bin_t = pdf_tbins-2;
  if (bin_t < 0)
    bin_t = 0;
  double t_d = (time_residual - (bin_t*this->pdf_deltat+this->pdf_mint))/(this->pdf_deltat);


  double pdf_val000 = this->pdf[(bin_s+0) * pdf_taubins * pdf_tbins + (bin_tau+0) * pdf_tbins + (bin_t+0)];
  double pdf_val100 = this->pdf[(bin_s+1) * pdf_taubins * pdf_tbins + (bin_tau+0) * pdf_tbins + (bin_t+0)];
  double pdf_val001 = this->pdf[(bin_s+0) * pdf_taubins * pdf_tbins + (bin_tau+0) * pdf_tbins + (bin_t+1)];
  double pdf_val101 = this->pdf[(bin_s+1) * pdf_taubins * pdf_tbins + (bin_tau+0) * pdf_tbins + (bin_t+1)];
  double pdf_val010 = this->pdf[(bin_s+0) * pdf_taubins * pdf_tbins + (bin_tau+1) * pdf_tbins + (bin_t+0)];
  double pdf_val110 = this->pdf[(bin_s+1) * pdf_taubins * pdf_tbins + (bin_tau+1) * pdf_tbins + (bin_t+0)];
  double pdf_val011 = this->pdf[(bin_s+0) * pdf_taubins * pdf_tbins + (bin_tau+1) * pdf_tbins + (bin_t+1)];
  double pdf_val111 = this->pdf[(bin_s+1) * pdf_taubins * pdf_tbins + (bin_tau+1) * pdf_tbins + (bin_t+1)];
  double pdf_val00 = pdf_val000*(1-s_d) + pdf_val100*s_d;
  double pdf_val01 = pdf_val001*(1-s_d) + pdf_val101*s_d;
  double pdf_val10 = pdf_val010*(1-s_d) + pdf_val110*s_d;
  double pdf_val11 = pdf_val011*(1-s_d) + pdf_val111*s_d;
  double pdf_val0 = pdf_val00*(1-tau_d) + pdf_val10*tau_d;
  double pdf_val1 = pdf_val01*(1-tau_d) + pdf_val11*tau_d;
  double pdf_val = pdf_val0*(1-t_d) + pdf_val1*t_d;

  return pdf_val;
}

