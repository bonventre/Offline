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
#include <TGraphAsymmErrors.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TCanvas.h>

#include "prototype/Utils.hh"
#include "fit/StrawFit.hh"

int main(int argc, char** argv)
{
  std::string filename = "";
  std::string outname = "";
  std::string help = "./efficiency_sim -f <input filename> -o <output filename> [-F (fix position)]";
  std::string inputcommand = "";
  int fix = 0;
  for (int i=0;i<argc;i++)
    inputcommand += std::string(argv[i]);
  
  int argc2 = 0;char **argv2;TApplication theApp("tapp", &argc2, argv2);

  int c;
  while ((c = getopt (argc, argv, "Fhf:o:")) != -1){
    switch (c){
      case 'F':
        fix = 1;
        break;
      case 'h':
        std::cout << help << std::endl;
        return 0;
      case 'f':
        filename = std::string(optarg);
        break;
      case 'o':
        outname = std::string(optarg);
        break;
      case '?':
        if (optopt == 'f' || optopt == 'o')
          std::cout << "Option -" << optopt << " requires an argument." << std::endl;
        else
          std::cout << "Unknown option `-" << (char)optopt << "'." << std::endl;
        return 1;
    }
  }

  if (filename.size() == 0 || outname.size() == 0){
    std::cout << help << std::endl;
    return 1;
  }
 
  gStyle->SetOptFit(1);

  std::vector<double> times;
  std::vector<double> tots;
  TFile *f = new TFile(filename.c_str());
  /*
  TDirectory *d = (TDirectory*) f->Get("SHD");
  TTree *t = (TTree*) d->Get("shdiag");
  Int_t straw, panel;
  Int_t mcproc;
  Float_t time[2], mcshd;
  Float_t ewm;
  Double_t mcsptime;
  Float_t tot[2];

  t->SetBranchAddress("straw",&straw);
  t->SetBranchAddress("panel",&panel);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("tot",&tot);
  t->SetBranchAddress("mcsptime",&mcsptime);
  t->SetBranchAddress("mcshd",&mcshd);
  t->SetBranchAddress("ewmoffset",&ewm);
  t->SetBranchAddress("mcproc",&mcproc);

  for (int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    if (mcproc != 56)
      continue;
    // the prototype is only the top 8 straws
    // and our generator is pointed at panel 1
    if (straw < 88 || panel != 1)
      continue;
    times.push_back(std::min(time[0],time[1]) - mcsptime + ewm);
    tots.push_back((tot[0]+tot[1])/2.);
  }
  */

  TDirectory *d = (TDirectory*) f->Get("makeSD");
  TTree *t = (TTree*) d->Get("sddiag");
  Int_t straw, panel;
  Int_t tdc[2], tot[2], mcproc;
  Double_t mctime;
  Float_t mcdca, ecptime, wdist[2], ewmoffset;
  t->SetBranchAddress("mcproc",&mcproc);
  t->SetBranchAddress("straw",&straw);
  t->SetBranchAddress("panel",&panel);
  t->SetBranchAddress("tdc",&tdc);
  t->SetBranchAddress("tot",&tot);
  t->SetBranchAddress("mctime",&mctime);
  t->SetBranchAddress("mcdca",&mcdca);
  t->SetBranchAddress("wdist",&wdist);
  t->SetBranchAddress("ecptime",&ecptime);
  t->SetBranchAddress("ewmoffset",&ewmoffset);


  for (int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    if (mcproc != 56)
      continue;
    if (straw < 88 || panel != 1)
      continue;
    times.push_back((tdc[0]+tdc[1])/2.0*2*0.015625 - (fmod(mctime,1695) + ecptime + ewmoffset));
    double totcal = tot[0] * 4 + (0x7F - tdc[0] & 0x7F)*0.015625*2;
    double tothv  = tot[1] * 4 + (0x7F - tdc[1] & 0x7F)*0.015625*2;
//    std::cout << "Cal: " << tot[0]*4 << " " << totcal << std::endl;
//    std::cout << "HV: " << tot[1]*4 << " " << tothv << std::endl;
    tots.push_back((totcal+tothv)/2.); 
  }
















  TH2F *h1 = new TH2F("h1","h1",32,0,64,120,-20,100);
  TH1F *h1d = new TH1F("h1d","h1d",32,0,64);
//  timeoffset = -20;

  for (int i=0;i<times.size();i++){
    h1->Fill(tots[i],times[i]);
  }

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TProfile *hp = (TProfile*) h1->ProfileX("hp",1,-1,"S");
  hp->SetLineColor(kBlue);
  hp->SetTitle("");
  hp->SetStats(0);
  hp->GetXaxis()->SetTitle("Time over threshold (ns)");
  hp->GetYaxis()->SetTitle("Drift time (ns)");
  hp->SetMarkerStyle(22);
  hp->GetXaxis()->SetRangeUser(0,50);
//  hp->Draw();
  h1->Draw();
//  h1d->Draw();

  theApp.Run();

  return 0;
}
