#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TError.h>
#include <TLine.h>
#include <TStyle.h>
#include <TFile.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TF1.h>

#include "prototype/Utils.hh"
#include "prototype/Parser.hh"
#include "event/Hit.h"

int main(int argc, char** argv)
{
  gErrorIgnoreLevel = 1001;
//  gROOT->SetBatch(kTRUE);

  int strawnum = -1;
  std::string filename_list = "";
  std::string outname = "";
  bool old_data = false;
  std::string help = "./longitudinal_data -s <strawnum> -f <filename list> -o <outname> [-u (old data)]";
  std::string inputcommand = "";
  for (int i=0;i<argc;i++)
    inputcommand += std::string(argv[i]);
  

  int c;
  while ((c = getopt (argc, argv, "hus:f:o:")) != -1){
    switch (c){
      case 'h':
        std::cout << help << std::endl;
        return 0;
      case 'f':
        filename_list = std::string(optarg);
        break;
      case 's':
        strawnum = atoi(optarg);
        break;
      case 'o':
        outname = std::string(optarg);
        break;
      case 'u':
        old_data = true;
        break;
      case '?':
        if (optopt == 'f' || optopt == 's' || optopt == 'o')
          std::cout << "Option -" << optopt << " requires an argument." << std::endl;
        else
          std::cout << "Unknown option `-" << optopt << "'." << std::endl;
        return 1;
    }
  }

  if (strawnum == -1 || filename_list.size() == 0 || outname.size() == 0){
    std::cout << help << std::endl;
    return 1;
  }

  std::vector<std::string> filenames;
  std::string line;
  std::ifstream fnlist(filename_list.c_str());
  while (std::getline(fnlist,line)){
    filenames.push_back(line);
  }

  int argc2 = 0;char **argv2;TApplication theApp("tapp", &argc2, argv2);
  
  gSystem->Load("/home/pixel/mu2e/newruns/software/Dict.so");

  gStyle->SetOptFit(1);

  std::vector<double> alltop_x; // all hits top pixel avg perpendicular to straw [pixel number (50um)]
  std::vector<double> allbot_x; // all hits bot pixel avg perpendicular to straw [pixel number (50um)]
  std::vector<double> alltop_y; // all hits top pixel avg along straw [pixel number (250um)]
  std::vector<double> allbot_y; // all hits bot pixel avg along straw [pixel number (250um)]


  std::vector<std::vector<bool> > has_straw(8, std::vector<bool>()); // has a hit in straw i
  std::vector<std::vector<bool> > has_shower(8, std::vector<bool>()); // has a hit in a straw not neighbors with straw i

  std::vector<std::vector<double> > alltimes(8, std::vector<double>()); // drift time measurement in straw i
  std::vector<std::vector<double> > alldts(8, std::vector<double>());
  std::vector<std::vector<double> > allmaxvals(8, std::vector<double>()); // NOT USED
  std::vector<std::vector<double> > alltots(8, std::vector<double>()); // NOT USED



  double total_time = 0; // total time length of these run before HV trips

  parse_data_files(filenames, strawnum, alltop_x, allbot_x, alltop_y, allbot_y, has_straw, has_shower, alltimes, alldts, allmaxvals, alltots, old_data);

  // make a vector with just the hits we want for this channel
  int istraw = strawnum;
  std::vector<double> dts;
  for (int j=0;j<alltop_y.size();j++){
    if (has_straw[istraw][j] && !has_shower[istraw][j]){
      dts.push_back(alldts[istraw][j]);
    }
  }

  // calculate mean so that it can be centered on zero
  std::vector<double> dts_copy = dts;
  std::sort(dts_copy.begin(),dts_copy.end());
  double median_dt = dts_copy[(int) (dts_copy.size()/2)];

  std::cout << "MEDIAN: " << median_dt << std::endl;

  TH1F *h = new TH1F("","",100,median_dt-10,median_dt+10);
  for (int i=0;i<dts.size();i++)
    h->Fill(dts[i]);
  TF1 *f1 = new TF1("f1","[2]*exp(-0.5*((x-[0])/[1])**2)",median_dt-10,median_dt+10);
  f1->SetParameter(2,h->GetMaximum());
  f1->SetParameter(0,median_dt);
  f1->SetParameter(1,h->GetRMS());
  h->Fit(f1,"QN0");



  double offset = f1->GetParameter(0);

  std::cout << "Mean: " << f1->GetParameter(0) << " +- " << f1->GetParError(0) << std::endl;
  std::cout << "Sigma: " << f1->GetParameter(1) << " +- " << f1->GetParError(1) << std::endl;

  //TCanvas *c1 = new TCanvas("c1","c1",800,600);
  //h->Draw();
  //f1->Draw("same");
 

  TH1F *h2 = new TH1F("","",300,-2,2);
  for (int i=0;i<dts.size();i++)
    h2->Fill(dts[i]-offset);
  f1->SetParameter(2,h2->GetMaximum());
  f1->SetParameter(0,0);
  f1->SetParameter(1,h2->GetRMS());
  h2->Fit(f1,"LQN0","",-1,1);

  offset += f1->GetParameter(0);
  f1->SetParameter(0,0);

  std::cout << "Mean: " << offset << " +- " << f1->GetParError(0) << std::endl;
  std::cout << "Sigma: " << f1->GetParameter(1) << " +- " << f1->GetParError(1) << std::endl;

  //TCanvas *c2 = new TCanvas("c2","c2",800,600);
  //h2->Draw();
  //f1->Draw("same");
  //theApp.Run();


 
  
  TH1F *h3 = new TH1F("","",75,-2,2);
  for (int i=0;i<dts.size();i++)
    h3->Fill(dts[i]-offset);

  h3->Draw();
  f1->SetParameter(0,0);
  f1->SetParameter(2,f1->GetParameter(2)*4);
  f1->SetRange(-1,1);
  f1->Draw("same");



  TFile *fout = new TFile(outname.c_str(),"RECREATE");
  TTree *tout = new TTree("longitudinal","longitudinal");
  double t_dt;
  double t_doca;
  tout->Branch("dt",&t_dt,"dt/D");
  for (int i=0;i<dts.size();i++){
    t_dt = dts[i];
    tout->Fill();
  }
  tout->Write();
  TTree *tresults = new TTree("longitudinalfit","longitudinalfit");
  double t_mean, t_sigma, t_scale;
  double t_meanerr, t_sigmaerr;
  tresults->Branch("mean",&t_mean,"mean/D");
  tresults->Branch("meanerr",&t_meanerr,"meanerr/D");
  tresults->Branch("sigma",&t_sigma,"sigma/D");
  tresults->Branch("sigmaerr",&t_sigmaerr,"sigmaerr/D");
  tresults->Branch("scale",&t_scale,"scale/D");
  t_mean = offset;
  t_sigma = f1->GetParameter(1);
  t_meanerr = f1->GetParError(0);
  t_sigmaerr = f1->GetParError(1);
  t_scale = f1->GetParameter(2);
  tresults->Fill();
  tresults->Write();

  TObjString *tl = new TObjString(inputcommand.c_str());
  fout->WriteObject(tl,"log");

  fout->Close();

  theApp.Run();
  
  return 0;
}
