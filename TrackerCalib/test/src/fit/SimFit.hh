#ifndef _SIMFIT_HH_
#define _SIMFIT_HH_

#include "fit/TrackerCalibFit.hh"

class SimFit : public TableFit {
  public:
    SimFit(std::vector<double> &_docas, std::vector<double> &_times,
        std::vector<double> &_constraint_means, std::vector<double> &_constraints,
        int _k, double voltage=1425.) : 
      TableFit(_k,_constraint_means, _constraints, voltage), docas(_docas), times(_times)
  {
    nparams = 4;
  };

    double operator() (const std::vector<double> &x) const;
    void calculate_weighted_pdf (const std::vector<double> &x, TH1F* h, double doca_min=-1, double doca_max=1e8, bool dist=false) const;

    std::vector<double> docas;
    std::vector<double> times;

};

#endif
