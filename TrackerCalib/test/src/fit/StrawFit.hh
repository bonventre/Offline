#ifndef _STRAWFIT_HH_
#define _STRAWFIT_HH_

#include <vector>
#include <TH1F.h>
#include <TFile.h>
#include <cmath>

#include <Minuit2/FCNBase.h>

#define STRAW_RADIUS 2.5
#define STRAW_INNER_RADIUS 2.485

class StrawFit : public ROOT::Minuit2::FCNBase {
  public:
    StrawFit(std::vector<double> &_alltops, std::vector<double> &_allbots, std::vector<bool> &_has_straw, std::vector<bool> &_has_shower) : alltops(_alltops), allbots(_allbots),has_straw(_has_straw), has_shower(_has_shower) {};

    double operator() (const std::vector<double> &x) const;
    double Up() const { return 0.5; };
    
    std::vector<double> alltops;
    std::vector<double> allbots;
    std::vector<bool> has_straw;
    std::vector<bool> has_shower;

};


class StrawFitPL : public ROOT::Minuit2::FCNBase {
  public:
    StrawFitPL(std::vector<double> &_alltops, std::vector<double> &_allbots, std::vector<bool> &_has_straw, std::vector<bool> &_has_shower) : alltops(_alltops), allbots(_allbots),has_straw(_has_straw), has_shower(_has_shower) {};

    double operator() (const std::vector<double> &x) const;
    double Up() const { return 0.5; };
    
    std::vector<double> alltops;
    std::vector<double> allbots;
    std::vector<bool> has_straw;
    std::vector<bool> has_shower;

};


class StrawFitPL2 : public ROOT::Minuit2::FCNBase {
  public:
    StrawFitPL2(std::vector<double> &_alltops, std::vector<double> &_allbots, std::vector<bool> &_has_straw, std::vector<bool> &_has_shower) : alltops(_alltops), allbots(_allbots),has_straw(_has_straw), has_shower(_has_shower) {
    
      TH1F* hpolya = new TH1F("hpolya","hpolya",500,0,10);
      hpcdf = new TH1F("hpcdf","hpcdf",500,0,10);

      for (int i=0;i<500;i++){
        double x = hpolya->GetBinCenter(i+1);
        double m = 0.25 + 1; // theta + 1
        double y = pow(m,m)/tgamma(m)*pow(x,m-1)*exp(-m*x);
        hpolya->SetBinContent(i+1,y);
        hpcdf->SetBinContent(i+1,hpcdf->GetBinContent(i)+hpolya->GetBinContent(i+1));
        hpcdf->SetBinContent(i+1,1-exp(-x));
      }
      hpcdf->Scale(1.0/hpcdf->GetBinContent(500));
    };

    double operator() (const std::vector<double> &x) const;
    double Up() const { return 0.5; };
    
    std::vector<double> alltops;
    std::vector<double> allbots;
    std::vector<bool> has_straw;
    std::vector<bool> has_shower;
    TH1F *hpcdf;

};


#endif
