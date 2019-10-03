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
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TApplication.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPaveText.h>

#include "fit/SimFit.hh"

#define LONG_VELOCITY 200

int main(int argc, char** argv)
{
  std::string filename[2];
  filename[0] = "";
  filename[1] = "";
  std::string help = "./compare_longitudinal -d <data filename> -s <sim filename>";
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

  if (filename[0].size() == 0 && filename[1].size() == 0 && optind == argc-2){
    for (int index = optind; index < argc; index++){
      filename[index-optind] = std::string(argv[index]);
    }
  }

  if (filename[0].size() == 0 || filename[1].size() == 0){
    std::cout << help << std::endl;
    return 1;
  }
 
  int argc2 = 0;char **argv2;TApplication theApp("tapp", &argc2, argv2);

  std::vector<double> dts[2];
  double mean[2], sigma[2], meanerr[2], sigmaerr[2], scale[2];

  TH1F* hresids[2];
  TF1* fpdfs[2];

  TH1F* hresidsd[2];
  TF1* fpdfsd[2];

  for (int ifile=0;ifile<2;ifile++){
    TFile *f = new TFile(filename[ifile].c_str());
    TTree *t = (TTree*) f->Get("longitudinal");
    double t_dt;
    t->SetBranchAddress("dt",&t_dt);
    for (int i=0;i<t->GetEntries();i++){
      t->GetEntry(i);
      dts[ifile].push_back(t_dt);
    }
    TTree *tfit = (TTree*) f->Get("longitudinalfit");
    tfit->SetBranchAddress("mean",&mean[ifile]);
    tfit->SetBranchAddress("sigma",&sigma[ifile]);
    tfit->SetBranchAddress("meanerr",&meanerr[ifile]);
    tfit->SetBranchAddress("sigmaerr",&sigmaerr[ifile]);
    tfit->SetBranchAddress("scale",&scale[ifile]);
    tfit->GetEntry(0);

    hresids[ifile] = new TH1F(TString::Format("hresid_%d",ifile),"hresid",75,-2,2);
    hresids[ifile]->Sumw2();
    hresidsd[ifile] = new TH1F(TString::Format("hresidd_%d",ifile),"hresid",75,-200,200);
    hresidsd[ifile]->Sumw2();
  
    fpdfs[ifile] = new TF1("","[2]*exp(-0.5*((x-[0])/[1])**2)",-1,1);
    fpdfs[ifile]->SetParameter(0,0);
    fpdfs[ifile]->SetParameter(1,sigma[ifile]);

    fpdfsd[ifile] = new TF1("","[2]*exp(-0.5*((x-[0])/[1])**2)",-100,100);
    fpdfsd[ifile]->SetParameter(0,0);
    fpdfsd[ifile]->SetParameter(1,sigma[ifile]*LONG_VELOCITY/2.);

    
    for (int i=0;i<dts[ifile].size();i++){
      hresids[ifile]->Fill(dts[ifile][i]-mean[ifile]);
      hresidsd[ifile]->Fill((dts[ifile][i]-mean[ifile])*LONG_VELOCITY/2.);
    }

    fpdfs[ifile]->SetParameter(2,scale[ifile]/hresids[ifile]->Integral());
    fpdfsd[ifile]->SetParameter(2,scale[ifile]/hresidsd[ifile]->Integral());
    
    hresids[ifile]->Scale(1.0/hresids[ifile]->Integral());
    hresidsd[ifile]->Scale(1.0/hresidsd[ifile]->Integral());
  }

  hresids[0]->SetStats(0);
  hresids[0]->SetTitle("");
  hresids[0]->GetXaxis()->SetTitle("Time Residual (ns)");
  hresidsd[0]->SetStats(0);
  hresidsd[0]->SetTitle("");
  hresidsd[0]->GetXaxis()->SetTitle("Longitudinal Distance Residual (mm)");


  fpdfs[0]->SetLineColor(kBlue);
  fpdfs[1]->SetLineColor(kRed);
  hresids[0]->SetLineColor(kBlue);
  hresids[1]->SetLineColor(kRed);

  TCanvas *cresid = new TCanvas("cresid","cresid",800,600);
  hresids[0]->Draw("E");
  hresids[1]->Draw("E same");
  fpdfs[0]->Draw("same");
  fpdfs[1]->Draw("same");

  TLegend *l = new TLegend(0.65,0.65,0.85,0.85);
  l->AddEntry(hresids[0],"Data","lepf");
  l->AddEntry(fpdfs[0],"Data fit","l");
  l->AddEntry(hresids[1],"Simulation","lepf");
  l->AddEntry(fpdfs[1],"Simulation fit","l");
  l->Draw();

  TPaveText *pt = new TPaveText(0.65,0.35,0.85,0.55,"NDC");
  TText *tt = pt->AddText(TString::Format("#sigma = %.3f #pm %.3f ns",sigma[0],sigmaerr[0]));
  tt->SetTextColor(kBlue);
  tt->Draw();
  tt = pt->AddText(TString::Format("#sigma = %.3f #pm %.3f ns",sigma[1],sigmaerr[1]));
  tt->SetTextColor(kRed);
  tt->Draw();
  pt->Draw();

  fpdfsd[0]->SetLineColor(kBlue);
  fpdfsd[1]->SetLineColor(kRed);
  hresidsd[0]->SetLineColor(kBlue);
  hresidsd[1]->SetLineColor(kRed);

  TCanvas *cresidd = new TCanvas("cresidd","cresidd",800,600);
  hresidsd[0]->Draw("E");
  hresidsd[1]->Draw("E same");
  fpdfsd[0]->Draw("same");
  fpdfsd[1]->Draw("same");

  TLegend *ld = new TLegend(0.65,0.65,0.85,0.85);
  ld->AddEntry(hresidsd[0],"Data","lepf");
  ld->AddEntry(fpdfsd[0],"Data fit","l");
  ld->AddEntry(hresidsd[1],"Simulation","lepf");
  ld->AddEntry(fpdfsd[1],"Simulation fit","l");
  ld->Draw();

  TPaveText *ptd = new TPaveText(0.65,0.35,0.85,0.55,"NDC");
  TText *ttd = ptd->AddText(TString::Format("#sigma = %.3f #pm %.3f mm",sigma[0]*LONG_VELOCITY/2.,sigmaerr[0]*LONG_VELOCITY/2.));
  ttd->SetTextColor(kBlue);
  ttd->Draw();
  ttd = ptd->AddText(TString::Format("#sigma = %.3f #pm %.3f mm",sigma[1]*LONG_VELOCITY/2.,sigmaerr[1]*LONG_VELOCITY/2.));
  ttd->SetTextColor(kRed);
  ttd->Draw();
  ptd->Draw();



  theApp.Run();

  return 0;

}

