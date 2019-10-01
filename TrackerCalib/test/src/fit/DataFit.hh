#ifndef _DATAFIT_HH_
#define _DATAFIT_HH_

#include "fit/TrackerCalibFit.hh"

class DataFit : public TableFit {
  public:
    DataFit(std::vector<double> &_tops, std::vector<double> &_bots, std::vector<double> &_times, std::vector<double> &_constraint_means, std::vector<double> &_constraints, int _k, double voltage=1425.)
      : TableFit(_k,_constraint_means,_constraints,voltage), tops(_tops), bots(_bots), times(_times)
    {
      nparams = 7;
    };

    double operator() (const std::vector<double> &x) const;
    
    void calculate_weighted_pdf (const std::vector<double> &x, TH1F* h, double doca_min=-1, double doca_max=1e8, bool dist=false) const;
    
    void metropolis(std::vector<double> &seed, std::vector<double> &errors, std::vector<std::vector<double> > &results, std::vector<double> &nlls, int steps);

    std::vector<double> tops;
    std::vector<double> bots;
    std::vector<double> times;
};

class GaussianDataFit : public Fit {
  public:
    GaussianDataFit(std::vector<double> &_tops, std::vector<double> &_bots, std::vector<double> &_times, std::vector<double> &_constraint_means, std::vector<double> &_constraints, double voltage=1425.) : Fit(_constraint_means,_constraints, voltage), tops(_tops), bots(_bots), times(_times){
      nparams = 7;
      outlier = std::vector<bool>(tops.size(),false);
    };

    double operator() (const std::vector<double> &x) const;
    
    void SelectOutliers(std::vector<double> &seed);

    std::vector<double> tops;
    std::vector<double> bots;
    std::vector<double> times;
    std::vector<bool> outlier;
};

#endif
