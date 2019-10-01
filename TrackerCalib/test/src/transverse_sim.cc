#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>

#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TMath.h>
#include <TObjString.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>

#include "fit/SimFit.hh"

int main(int argc, char** argv)
{
  float voltage = 1400;
  std::string filename = "";
  std::string outname = "";
  std::string help = "./transverse_sim -f <filename list> -o <outname> [-v voltage]";
  std::string inputcommand = "";
  for (int i=0;i<argc;i++)
    inputcommand += std::string(argv[i]);
  

  int c;
  while ((c = getopt (argc, argv, "hf:o:v:")) != -1){
    switch (c){
      case 'h':
        std::cout << help << std::endl;
        return 0;
      case 'v':
        voltage = atof(optarg);      
        break;
      case 'f':
        filename = std::string(optarg);
        break;
      case 'o':
        outname = std::string(optarg);
        break;
      case '?':
        if (optopt == 'v' || optopt == 'f' || optopt == 'o')
          std::cout << "Option -" << optopt << " requires an argument." << std::endl;
        else
          std::cout << "Unknown option `-" << optopt << "'." << std::endl;
        return 1;
    }
  }

  if (filename.size() == 0 || outname.size() == 0){
    std::cout << help << std::endl;
    return 1;
  }
 
  TFile *f = new TFile(filename.c_str());
  TDirectory *d = (TDirectory*) f->Get("SHD");
  TTree *t = (TTree*) d->Get("shdiag");
  Int_t straw, panel;
  Int_t mcproc;
  Float_t time[2], mcshd;
  Float_t ewm;
  Double_t mcsptime;

  t->SetBranchAddress("straw",&straw);
  t->SetBranchAddress("panel",&panel);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("mcsptime",&mcsptime);
  t->SetBranchAddress("mcshd",&mcshd);
  t->SetBranchAddress("ewmoffset",&ewm);
  t->SetBranchAddress("mcproc",&mcproc);

  std::vector<double> times;
  std::vector<double> docas;
  for (int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    if (mcproc != 56)
      continue;
    // the prototype is only the top 8 straws
    // and our generator is pointed at panel 1
    if (straw < 88 || panel != 1)
      continue;
    times.push_back(std::min(time[0],time[1]) - mcsptime + ewm);
    docas.push_back(fabs(mcshd));
  }

  /*
  TDirectory *d = (TDirectory*) f->Get("makeSD");
  TTree *t = (TTree*) d->Get("sddiag");
  Int_t straw, panel;
  Int_t mcproc;
  Int_t tdc[2];
  Double_t mctime;
  Float_t mcdca, ecptime, wdist[2];
  Float_t ewm;

  t->SetBranchAddress("straw",&straw);
  t->SetBranchAddress("panel",&panel);
  t->SetBranchAddress("tdc",&tdc);
  t->SetBranchAddress("mctime",&mctime);
  t->SetBranchAddress("mcdca",&mcdca);
  t->SetBranchAddress("wdist",&wdist);
  t->SetBranchAddress("ewmoffset",&ewm);
  t->SetBranchAddress("mcproc",&mcproc);

  std::vector<double> times;
  std::vector<double> docas;
  for (int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    if (mcproc != 56)
      continue;
    // the prototype is only the top 8 straws
    // and our generator is pointed at panel 1
    if (straw < 88 || panel != 1)
      continue;
    times.push_back(std::min(tdc[0],tdc[1])*2*0.015625 - (fmod(mctime,1695) + ewm - 10));
    docas.push_back(fabs(mcdca));
  }
  */



  std::cout << times.size() << " total events" << std::endl;

  // x[0] is the tau
  // x[1] is the gaussian width
  // x[2] is time offset between PMT time and straw hit time, so time[i] = x[2] when the ion is directly on the wire
  // x[3] is the background
  std::vector<double> seed(4,0);
  std::vector<double> errors(4,0);
  seed[0] = 10;
  seed[1] = 1;
  seed[2] = 0;
  seed[3] = 0;
  errors[0] = 0.1;
  errors[1] = 0.1;
  errors[2] = 0.1;
  errors[3] = 0.1;

  std::vector<double> constraint_means(4,0);
  std::vector<double> constraints(4,0);

  SimFit fit(docas,times,constraint_means,constraints,1,voltage);

  ROOT::Minuit2::MnStrategy mnStrategy(2);
  ROOT::Minuit2::MnUserParameters params(seed,errors);
  ROOT::Minuit2::MnMigrad migrad(fit,params,mnStrategy);

  migrad.SetLimits((unsigned) 0, fit.pdf_mintau,fit.pdf_maxtau);
  migrad.SetLimits((unsigned) 1, fit.pdf_mins,fit.pdf_maxs);
  migrad.Fix((unsigned) 3);

  ROOT::Minuit2::FunctionMinimum min = migrad();
  ROOT::Minuit2::MnUserParameters results = min.UserParameters();
  double minval = min.Fval();
  std::vector<double> bestfit = results.Params();
  std::vector<double> bestfiterrors = results.Errors();

  std::vector<std::string> names;
  names.push_back("tau (ns)");
  names.push_back("sigma (ns)");
  names.push_back("time offset (ns)");
  names.push_back("background (frac)");

  std::cout << "NLL: " << minval << std::endl;
  for (int i=0;i<names.size();i++){
    std::cout << i << names[i] << " : " << bestfit[i] << " +- " << bestfiterrors[i] << std::endl;
  }

  TFile *fout = new TFile(outname.c_str(),"RECREATE");
  TTree *tout = new TTree("transverse","transverse");
  double t_time;
  double t_doca;
  tout->Branch("time",&t_time,"time/D");
  tout->Branch("doca",&t_doca,"doca/D");
  for (int i=0;i<times.size();i++){
    t_time = times[i];
    t_doca = docas[i];
    tout->Fill();
  }
  tout->Write();
  TTree *tresults = new TTree("transversefit","transversefit");
  double t_tau, t_sigma, t_timeoffset, t_background;
  double t_tauerr, t_sigmaerr, t_timeoffseterr, t_backgrounderr;
  double t_nll, t_voltage;
  tresults->Branch("nll",&t_nll,"nll/D");
  tresults->Branch("voltage",&t_voltage,"voltage/D");
  tresults->Branch("tau",&t_tau,"tau/D");
  tresults->Branch("sigma",&t_sigma,"sigma/D");
  tresults->Branch("timeoffset",&t_timeoffset,"timeoffset/D");
  tresults->Branch("background",&t_background,"background/D");
  tresults->Branch("tauerr",&t_tauerr,"tauerr/D");
  tresults->Branch("sigmaerr",&t_sigmaerr,"sigmaerr/D");
  tresults->Branch("timeoffseterr",&t_timeoffseterr,"timeoffseterr/D");
  tresults->Branch("backgrounderr",&t_backgrounderr,"backgrounderr/D");
  t_nll = minval;
  t_voltage = voltage;
  t_tau = bestfit[0];
  t_sigma = bestfit[1];
  t_timeoffset = bestfit[2];
  t_background = bestfit[3];
  t_tauerr = bestfiterrors[0];
  t_sigmaerr = bestfiterrors[1];
  t_timeoffseterr = bestfiterrors[2];
  t_backgrounderr = bestfiterrors[3];
  tresults->Fill();
  tresults->Write();

  TObjString *tl = new TObjString(inputcommand.c_str());
  fout->WriteObject(tl,"log");
  fout->Close();


  return 0;
}
