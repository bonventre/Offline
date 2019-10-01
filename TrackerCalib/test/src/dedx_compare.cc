#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TApplication.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TF1.h>

#include "fit/StrawFit.hh"

int main(int argc, char** argv)
{
  std::string filename[2];
  filename[0] = "";
  filename[1] = "";
  std::string help = "./compare_dedx -d <data efficiency fit> -s <sim efficiency fit>";
  std::string inputcommand = "";
  for (int i=0;i<argc;i++)
    inputcommand += std::string(argv[i]);
  

  int c;
  while ((c = getopt (argc, argv, "hd:s:")) != -1){
    switch (c){
      case 'h':
        std::cout << help << std::endl;
        return 0;
      case 'd':
        filename[0] = std::string(optarg);      
        break;
      case 's':
        filename[1] = std::string(optarg);
        break;
      case '?':
        if (optopt == 'd' || optopt == 's')
          std::cout << "Option -" << optopt << " requires an argument." << std::endl;
        else
          std::cout << "Unknown option `-" << optopt << "'." << std::endl;
        return 1;
    }
  }

  if (filename[0].size() == 0 && filename[1].size() == 0 && optind-argc <= 2){
    for (int index = optind; index < argc; index++){
      filename[index-optind] = std::string(argv[index]);
    }
  }

  if (filename[0].size() == 0 && filename[1].size() == 0){
    std::cout << help << std::endl;
    return 1;
  }

  int filecount = 1;
  if (filename[1].size() > 0)
    filecount = 2;

  int argc2 = 0;char **argv2;TApplication theApp("tapp", &argc2, argv2);

  std::vector<double> docas[2];
  std::vector<double> pmps[2];

  TH2F *hdedx[2];
  TH1F *hprofile[2];
  TH1F *hprofileS[2];

  for (int ifile=0;ifile<filecount;ifile++){
    TFile *f = new TFile(filename[ifile].c_str());
    TTree *t = (TTree*) f->Get("efficiency");
    double t_doca;
    double t_pmp;
    t->SetBranchAddress("doca",&t_doca);
    t->SetBranchAddress("pmp",&t_pmp);
    for (int i=0;i<t->GetEntries();i++){
      t->GetEntry(i);
      docas[ifile].push_back(t_doca);
      pmps[ifile].push_back(t_pmp);
    }
    
    hdedx[ifile] = new TH2F(TString::Format("hdedx_%d",ifile),"hdedx",120,-1,5,150,0,1500);
    for (int j=0;j<docas[ifile].size();j++){
      double pathlength;
      if (docas[ifile][j] < STRAW_INNER_RADIUS)
        pathlength = sqrt(STRAW_INNER_RADIUS*STRAW_INNER_RADIUS-docas[ifile][j]*docas[ifile][j])*2;
      else
        pathlength = STRAW_INNER_RADIUS-docas[ifile][j];
      hdedx[ifile]->Fill(pathlength,pmps[ifile][j]*1500/4096.);
    }
    hprofile[ifile] = (TH1F*) hdedx[ifile]->ProfileX(TString::Format("hprofile_%d",ifile),-1,999);
    hprofileS[ifile] = (TH1F*) hdedx[ifile]->ProfileX(TString::Format("hprofileS_%d",ifile),-1,999,"S");
  }

  hprofile[0]->SetTitle("");
  hprofile[0]->SetStats(0);
  hprofile[0]->SetLineColor(kBlue);
  hprofile[0]->GetXaxis()->SetTitle("Path length through straw (mm)");
  hprofile[0]->GetYaxis()->SetTitle("Average peak minus pedestal (V)");
  hprofile[0]->GetXaxis()->SetRangeUser(0,5);
  hprofile[0]->GetYaxis()->SetRangeUser(0,250);

  hprofile[0]->Fit("pol1");

  double intercept = hprofile[0]->GetFunction("pol1")->GetParameter(0);
  double slope = hprofile[0]->GetFunction("pol1")->GetParameter(1);
  double intercepterr = hprofile[0]->GetFunction("pol1")->GetParError(0);
  double slopeerr = hprofile[0]->GetFunction("pol1")->GetParError(1);

  std::cout << intercept << " +- " << intercepterr << std::endl; 
  std::cout << slope << " +- " << slopeerr << std::endl; 

  TCanvas* cdedx = new TCanvas("cdedx","cdedx",800,600);
  hprofile[0]->Draw("E");
  TLegend *l = new TLegend(0.65,0.65,0.85,0.85);
  if (filecount > 1){
    hprofile[1]->SetLineColor(kRed);
    hprofile[1]->Draw("E same");

    hprofile[1]->Fit("pol1");
    
    double intercept1 = hprofile[1]->GetFunction("pol1")->GetParameter(0);
    double slope1 = hprofile[1]->GetFunction("pol1")->GetParameter(1);
    double intercepterr1 = hprofile[1]->GetFunction("pol1")->GetParError(0);
    double slopeerr1 = hprofile[1]->GetFunction("pol1")->GetParError(1);

    std::cout << intercept1 << " +- " << intercepterr1 << std::endl; 
    std::cout << slope1 << " +- " << slopeerr1 << std::endl; 



    l->AddEntry(hprofile[0],"Data","lepf");
    l->AddEntry(hprofile[1],"Simulation","lepf");
    l->Draw();
  }

  hprofileS[0]->SetTitle("");
  hprofileS[0]->SetStats(0);
  hprofileS[0]->SetLineColor(kBlue);
  hprofileS[0]->GetXaxis()->SetTitle("Path length through straw (mm)");
  hprofileS[0]->GetYaxis()->SetTitle("Average peak minus pedestal (V)");
  hprofileS[0]->GetXaxis()->SetRangeUser(0,5);
  hprofileS[0]->GetYaxis()->SetRangeUser(0,1000);

  TCanvas* cdedxS = new TCanvas("cdedxS","cdedx",800,600);
  hprofileS[0]->Draw("E");
  TLegend *lS = new TLegend(0.65,0.65,0.85,0.85);
  if (filecount > 1){
    hprofileS[1]->SetLineColor(kRed);
    hprofileS[1]->Draw("E same");

    lS->AddEntry(hprofileS[0],"Data","lepf");
    lS->AddEntry(hprofileS[1],"Simulation","lepf");
    lS->Draw();
  }

  TH1F *hdedxerr[2];
  hdedxerr[0] = new TH1F(TString::Format("hdedxerr_%d",0),"hdedx",120,-1,5);
  hdedxerr[1] = new TH1F(TString::Format("hdedxerr_%d",1),"hdedx",120,-1,5);
  for (int i=0;i<120;i++){
    hdedxerr[0]->SetBinContent(i+1,hprofileS[0]->GetBinError(i+1));
    if (filecount > 1)
      hdedxerr[1]->SetBinContent(i+1,hprofileS[1]->GetBinError(i+1));
  }

  TCanvas* cdedxerr = new TCanvas("cdedxerr","cdedxerr",800,600);
  hdedxerr[0]->SetTitle("");
  hdedxerr[0]->SetStats(0);
  hdedxerr[0]->SetLineColor(kBlue);
  hdedxerr[0]->GetXaxis()->SetTitle("Path length through straw (mm)");
  hdedxerr[0]->GetYaxis()->SetTitle("RMS of peak minus pedestal (V)");
  hdedxerr[0]->GetXaxis()->SetRangeUser(0,5);
  hdedxerr[0]->GetYaxis()->SetRangeUser(0,1000);

  hdedxerr[0]->Draw("");
  TLegend *lerr = new TLegend(0.65,0.65,0.85,0.85);
   if (filecount > 1){
    hdedxerr[1]->SetLineColor(kRed);
    hdedxerr[1]->Draw("same");

    lerr->AddEntry(hdedxerr[0],"Data","l");
    lerr->AddEntry(hdedxerr[1],"Simulation","l");
    lerr->Draw();
  }






  theApp.Run();

  return 0;

}
