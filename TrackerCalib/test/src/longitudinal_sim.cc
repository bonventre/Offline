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

int main(int argc, char** argv)
{
  std::string filename = "";
  std::string outname = "";
  std::string help = "./longitudinal_sim -f <filename list> -o <outname>";
  std::string inputcommand = "";
  for (int i=0;i<argc;i++)
    inputcommand += std::string(argv[i]);
  

  int c;
  while ((c = getopt (argc, argv, "hf:o:")) != -1){
    switch (c){
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
          std::cout << "Unknown option `-" << optopt << "'." << std::endl;
        return 1;
    }
  }

  if (filename.size() == 0 || outname.size() == 0){
    std::cout << help << std::endl;
    return 1;
  }
 
  int argc2 = 0;char **argv2;TApplication theApp("tapp", &argc2, argv2);

  gStyle->SetOptFit(1);

  TFile *f = new TFile(filename.c_str());
  TDirectory *d = (TDirectory*) f->Get("SHD");
  TTree *t = (TTree*) d->Get("shdiag");
  Float_t time[2];
  Float_t mcshlen;
  Int_t straw, panel, mcproc;

  t->SetBranchAddress("straw",&straw);
  t->SetBranchAddress("panel",&panel);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("mcproc",&mcproc);
  t->SetBranchAddress("mcshlen",&mcshlen);

  std::vector<double> dts;
  std::vector<double> dists;
  for (int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    if (mcproc != 56)
      continue;
    // the prototype is only the top 8 straws
    // and our generator is pointed at panel 1
    if (straw < 88 || panel != 1)
      continue;
    dts.push_back(time[0]-time[1]);
    dists.push_back(-1*mcshlen);
  }

  std::cout << dts.size() << " total events" << std::endl;



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

  TH1F *h2 = new TH1F("","",300,-2,2);
  for (int i=0;i<dts.size();i++)
    h2->Fill(dts[i]-offset);
  f1->SetParameter(0,0);
  f1->SetParameter(2,h2->GetMaximum());
  f1->SetParameter(1,h2->GetRMS());
  h2->Fit(f1,"LQN0","",-1,1);

  offset += f1->GetParameter(0);
  f1->SetParameter(0,0);

  std::cout << "Mean: " << offset << " +- " << f1->GetParError(0) << std::endl;
  std::cout << "Sigma: " << f1->GetParameter(1) << " +- " << f1->GetParError(1) << std::endl;

  
  TH1F *h3 = new TH1F("","",75,-2,2);
  for (int i=0;i<dts.size();i++)
    h3->Fill(dts[i]-offset);

  h3->Draw();
  f1->SetParameter(0,0);
  f1->SetParameter(2,f1->GetParameter(2)*300/75.);
  //f1->SetParameter(2,h3->GetMaximum());
  //f1->SetParameter(1,h3->GetRMS());
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
