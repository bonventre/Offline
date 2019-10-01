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

#include "fit/SimFit.hh"

int main(int argc, char** argv)
{
  std::string filename[2];
  filename[0] = "";
  filename[1] = "";
  std::string help = "./compare_transverse -d <data filename> -s <sim filename>";
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

  std::vector<double> times[2];
  std::vector<double> docas[2];
  double tau[2], sigma[2], timeoffset[2], background[2];
  double tauerr[2], sigmaerr[2], timeoffseterr[2], backgrounderr[2];
  double voltage[2];
  std::vector<double> bestfit[2];
  SimFit *fits[2];
  std::vector<double> constraint_means(4,0);
  std::vector<double> constraints(4,0);

  TH1F* hresids[2];
  TH1F* hpdfs[2];

  TH1F* hresidsd[2];
  TH1F* hpdfsd[2];

  TH1F* weights[2];

  for (int ifile=0;ifile<2;ifile++){
    TFile *f = new TFile(filename[ifile].c_str());
    TTree *t = (TTree*) f->Get("transverse");
    double t_time;
    double t_doca;
    t->SetBranchAddress("time",&t_time);
    t->SetBranchAddress("doca",&t_doca);
    for (int i=0;i<t->GetEntries();i++){
      t->GetEntry(i);
      times[ifile].push_back(t_time);
      docas[ifile].push_back(t_doca);
    }
    TTree *tfit = (TTree*) f->Get("transversefit");
    tfit->SetBranchAddress("voltage",&voltage[ifile]);
    tfit->SetBranchAddress("tau",&tau[ifile]);
    tfit->SetBranchAddress("sigma",&sigma[ifile]);
    tfit->SetBranchAddress("timeoffset",&timeoffset[ifile]);
    tfit->SetBranchAddress("background",&background[ifile]);
    tfit->SetBranchAddress("tauerr",&tauerr[ifile]);
    tfit->SetBranchAddress("sigmaerr",&sigmaerr[ifile]);
    tfit->SetBranchAddress("timeoffseterr",&timeoffseterr[ifile]);
    tfit->SetBranchAddress("backgrounderr",&backgrounderr[ifile]);
    tfit->GetEntry(0);

    bestfit[ifile] = std::vector<double>(4,0);
    bestfit[ifile][0] = tau[ifile];
    bestfit[ifile][1] = sigma[ifile];
    bestfit[ifile][2] = timeoffset[ifile];
    bestfit[ifile][3] = background[ifile];

    fits[ifile] = new SimFit(docas[ifile],times[ifile],constraint_means,constraints,1,voltage[ifile]);

    hresids[ifile] = new TH1F(TString::Format("hresid_%d",ifile),"hresid",200,-20,20);
    hpdfs[ifile] = new TH1F(TString::Format("hpdf_%d",ifile),"hpdf",200,-20,20);
    hresids[ifile]->Sumw2();
    hresidsd[ifile] = new TH1F(TString::Format("hresidd_%d",ifile),"hresid",125,-0.5,2);
    hpdfsd[ifile] = new TH1F(TString::Format("hpdfd_%d",ifile),"hpdf",125,-0.5,2);
    hresidsd[ifile]->Sumw2();

    weights[ifile] = new TH1F(TString::Format("weights_%d",ifile),"weights",10,0,2.5);
    for (int i=0;i<docas[ifile].size();i++){
      if (docas[ifile][i] > 2.5)
        continue;
      weights[ifile]->Fill(docas[ifile][i]);
    }
  }

  if (voltage[0] != voltage[1]){
    std::cout << "WARNING: looks like data and sim were fit at different HV settings: " << voltage[0] << " " << voltage[1] << std::endl;
  }

  // we want the distribution in doca to be the same 
  // in data and simulation, so we will accept a subset
  // of mc events to enforce this
  weights[1]->Divide(weights[0]);
  double mc_scale = weights[1]->GetMinimum();
  weights[1]->Scale(0);
  weights[0]->Scale(mc_scale);

  std::vector<double> times_temp = times[1];
  std::vector<double> docas_temp = docas[1];
  times[1].clear();
  docas[1].clear();
  for (int i=0;i<times_temp.size();i++){
    int bin = weights[1]->FindBin(docas_temp[i]);
    if (weights[1]->GetBinContent(bin) < weights[0]->GetBinContent(bin)){
      times[1].push_back(times_temp[i]);
      docas[1].push_back(docas_temp[i]);
      weights[1]->Fill(docas_temp[i]);
    }
  }

  for (int ifile=0;ifile<2;ifile++){

    fits[ifile]->calculate_weighted_pdf(bestfit[ifile],hpdfs[ifile]);
    fits[ifile]->calculate_weighted_pdf(bestfit[ifile],hpdfsd[ifile],-1,1e8,true);

    for (int i=0;i<docas[ifile].size();i++){
      double tr = fits[ifile]->TimeResidual(docas[ifile][i], times[ifile][i], timeoffset[ifile]);
      hresids[ifile]->Fill(tr);
      double doca_reco = fits[ifile]->T2D(times[ifile][i], timeoffset[ifile]);
      hresidsd[ifile]->Fill(doca_reco-docas[ifile][i]);
    }
    hresids[ifile]->Scale(1.0/hresids[ifile]->Integral());
    hpdfs[ifile]->Scale(hresids[ifile]->Integral()/hpdfs[ifile]->Integral());
    hresidsd[ifile]->Scale(1.0/hresidsd[ifile]->Integral());
    hpdfsd[ifile]->Scale(hresidsd[ifile]->Integral()/hpdfsd[ifile]->Integral());
  }

  hresids[0]->SetStats(0);
  hresids[0]->SetTitle("");
  hresids[0]->GetXaxis()->SetTitle("Time Residual (ns)");
  hresidsd[0]->SetStats(0);
  hresidsd[0]->SetTitle("");
  hresidsd[0]->GetXaxis()->SetTitle("Drift Radius Residual (mm)");


  hpdfs[0]->SetLineColor(kBlue);
  hpdfs[1]->SetLineColor(kRed);
  hresids[0]->SetLineColor(kBlue);
  hresids[1]->SetLineColor(kRed);

  TCanvas *cresid = new TCanvas("cresid","cresid",800,600);
  hresids[0]->Draw("E");
  hresids[1]->Draw("E same");
  hpdfs[0]->Draw("same");
  hpdfs[1]->Draw("same");

  TLegend *l = new TLegend(0.65,0.65,0.85,0.85);
  l->AddEntry(hresids[0],"Data","lepf");
  l->AddEntry(hpdfs[0],"Data fit","l");
  l->AddEntry(hresids[1],"Simulation","lepf");
  l->AddEntry(hpdfs[1],"Simulation fit","l");
  l->Draw();

  TPaveText *pt = new TPaveText(0.65,0.35,0.85,0.55,"NDC");
  TText *tt = pt->AddText(TString::Format("#sigma = %.3f #pm %.3f ns",sigma[0],sigmaerr[0]));
  tt->SetTextColor(kBlue);
  tt->Draw();
  tt = pt->AddText(TString::Format("#tau = %.3f #pm %.3f ns",tau[0],tauerr[0]));
  tt->SetTextColor(kBlue);
  tt->Draw();
  tt = pt->AddText(TString::Format("#sigma = %.3f #pm %.3f ns",sigma[1],sigmaerr[1]));
  tt->SetTextColor(kRed);
  tt->Draw();
  tt = pt->AddText(TString::Format("#tau = %.3f #pm %.3f ns",tau[1],tauerr[1]));
  tt->SetTextColor(kRed);
  tt->Draw();
  pt->Draw();

  hpdfsd[0]->SetLineColor(kBlue);
  hpdfsd[1]->SetLineColor(kRed);
  hresidsd[0]->SetLineColor(kBlue);
  hresidsd[1]->SetLineColor(kRed);

  TCanvas *cresidd = new TCanvas("cresidd","cresidd",800,600);
  hresidsd[0]->Draw("E");
  hresidsd[1]->Draw("E same");
  hpdfsd[0]->Draw("same");
  hpdfsd[1]->Draw("same");

  TLegend *ld = new TLegend(0.65,0.65,0.85,0.85);
  ld->AddEntry(hresidsd[0],"Data","lepf");
  ld->AddEntry(hpdfsd[0],"Data fit","l");
  ld->AddEntry(hresidsd[1],"Simulation","lepf");
  ld->AddEntry(hpdfsd[1],"Simulation fit","l");
  ld->Draw();

  TPaveText *ptd = new TPaveText(0.65,0.35,0.85,0.55,"NDC");
  TText *ttd = ptd->AddText(TString::Format("#sigma = %.3f #pm %.3f mm",sigma[0]*AVG_VELOCITY,sigmaerr[0]*AVG_VELOCITY));
  ttd->SetTextColor(kBlue);
  ttd->Draw();
  ttd = ptd->AddText(TString::Format("#tau = %.3f #pm %.3f mm",tau[0]*AVG_VELOCITY,tauerr[0]*AVG_VELOCITY));
  ttd->SetTextColor(kBlue);
  ttd->Draw();
  ttd = ptd->AddText(TString::Format("#sigma = %.3f #pm %.3f mm",sigma[1]*AVG_VELOCITY,sigmaerr[1]*AVG_VELOCITY));
  ttd->SetTextColor(kRed);
  ttd->Draw();
  ttd = ptd->AddText(TString::Format("#tau = %.3f #pm %.3f mm",tau[1]*AVG_VELOCITY,tauerr[1]*AVG_VELOCITY));
  ttd->SetTextColor(kRed);
  ttd->Draw();
  ptd->Draw();



  theApp.Run();

  return 0;

}

