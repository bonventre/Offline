// Base classes for fitting and for tabular drift PDF

#ifndef _STRAWCALIBFIT_HH_
#define _STRAWCALIBFIT_HH_

#include <vector>

#include <Minuit2/FCNBase.h>
#include <TH1F.h>

#include "fit/StrawDrift.hh"

#define pdf_tbins 500
#define pdf_taubins 50
#define pdf_sbins 25

#define AVG_VELOCITY 0.065

class Fit : public ROOT::Minuit2::FCNBase {
  public:
    Fit(std::vector<double> &_constraint_means, std::vector<double> &_constraints, double voltage=1425.) : constraint_means(_constraint_means), constraints(_constraints) {
      sd = new StrawDrift("src/fit/E2v.tbl",voltage,3000);
      nparams = 0;
    };

    double operator() (const std::vector<double> &x) const { return 0; };
    double Up() const { return 0.5; };

    double D2T(double distance, double time_offset) const;
    double T2D(double time, double time_offset) const;
    double TimeResidual(double doca, double time, double time_offset) const;
    
    std::vector<double> constraint_means;
    std::vector<double> constraints;
    StrawDrift *sd;
    int nparams;

};

class TableFit : public Fit {
  public:
    TableFit(int _k, std::vector<double> &_constraint_means, std::vector<double> &_constraints, double voltage=1425.);

    void calculate_full_pdf();

    double interpolatePDF(double time_residual, double sigma, double tau) const;

    double pdf_mint, pdf_mintau, pdf_mins;
    double pdf_maxt, pdf_maxtau, pdf_maxs;
    double pdf_deltat, pdf_deltatau, pdf_deltas;
    double *pdf_sigmas, *pdf_taus, *pdf_times;
    double *pdf;
    int k;
};

#endif

