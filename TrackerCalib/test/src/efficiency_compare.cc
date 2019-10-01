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
#include <TMath.h>

#include "prototype/Utils.hh"
#include "fit/StrawFit.hh"

int main(int argc, char** argv)
{
  std::string filename[2];
  bool recalculate = false;
  filename[0] = "";
  filename[1] = "";
  std::string help = "./compare_efficiency -d <data filename> -s <sim filename> [-r recalculate]";
  std::string inputcommand = "";
  for (int i=0;i<argc;i++)
    inputcommand += std::string(argv[i]);
  

  int c;
  while ((c = getopt (argc, argv, "hrd:s:")) != -1){
    switch (c){
      case 'h':
        std::cout << help << std::endl;
        return 0;
      case 'r':
        recalculate = true;
        break;
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

  std::vector<double> top[2];
  std::vector<double> bot[2];
  std::vector<double> doca[2];
  std::vector<bool> has_straw[2];
  std::vector<bool> has_shower[2];

  double meaneff[2], exponent[2], hdist[2], vdist[2];
  double meanefferr[2], exponenterr[2], hdisterr[2], vdisterr[2];

  TH1F* hefficiency[2];
  TH1F* hefficiency_d[2];
  TH1F* hpdf[2];
  TLine* lavg[2];

  double effavg[2], effavgerr[2];

  for (int ifile=0;ifile<2;ifile++){
    TFile *f = new TFile(filename[ifile].c_str());
    TTree *t = (TTree*) f->Get("efficiency");
    double t_top, t_bot, t_doca;
    int t_hit, t_shower;
    t->SetBranchAddress("top",&t_top);
    t->SetBranchAddress("bot",&t_bot);
    t->SetBranchAddress("doca",&t_doca);
    t->SetBranchAddress("hit",&t_hit);
    t->SetBranchAddress("shower",&t_shower);
    for (int i=0;i<t->GetEntries();i++){
      t->GetEntry(i);
      top[ifile].push_back(t_top);
      bot[ifile].push_back(t_bot);
      doca[ifile].push_back(t_doca);
      has_straw[ifile].push_back((t_hit == 1));
      has_shower[ifile].push_back((t_shower == 1));
    }

    TTree *tfit = (TTree*) f->Get("efficiencyfit");
    tfit->SetBranchAddress("meaneff",&meaneff[ifile]);
    tfit->SetBranchAddress("exponent",&exponent[ifile]);
    tfit->SetBranchAddress("hdist",&hdist[ifile]);
    tfit->SetBranchAddress("vdist",&vdist[ifile]);
    tfit->SetBranchAddress("meanefferr",&meanefferr[ifile]);
    tfit->SetBranchAddress("exponenterr",&exponenterr[ifile]);
    tfit->SetBranchAddress("hdisterr",&hdisterr[ifile]);
    tfit->SetBranchAddress("vdisterr",&vdisterr[ifile]);
    tfit->GetEntry(0);


    hefficiency[ifile] = new TH1F(TString::Format("eff_%d",ifile),"eff",80,0,4);
    hefficiency_d[ifile] = new TH1F(TString::Format("eff_d_%d",ifile),"eff",80,0,4);

    for (int j=0;j<top[ifile].size();j++){
      if (!has_shower[ifile][j]){
        double straw_dist;
        if (recalculate)
          straw_dist = calculate_DOCA(top[ifile][j],bot[ifile][j],hdist[ifile],vdist[ifile]);
        else
          straw_dist = doca[ifile][j];
        hefficiency_d[ifile]->Fill(straw_dist);
        if (has_straw[ifile][j]){
          hefficiency[ifile]->Fill(straw_dist);
        }
      }
    }

    hpdf[ifile] = new TH1F(TString::Format("pdf_%d",ifile),"pdf",80,0,4);
    for (int j=0;j<hpdf[ifile]->GetNbinsX();j++){
      double dist = hpdf[ifile]->GetBinCenter(j+1);
      double pdf_prob = meaneff[ifile] * (1 - 0.5*(1+TMath::Erf((dist-2.5)/(fabs(exponent[ifile])*sqrt(2)))));
      hpdf[ifile]->SetBinContent(j+1,pdf_prob);
    }

    double tot_d = hefficiency_d[ifile]->Integral(1,40);
    double tot_n = hefficiency[ifile]->Integral(1,40);
    effavg[ifile] = tot_n/tot_d;
    effavgerr[ifile] = sqrt(effavg[ifile]*(1-effavg[ifile])/tot_d);

    lavg[ifile] = new TLine(0,effavg[ifile],2.0,effavg[ifile]);

    for (int i=0;i<hefficiency[ifile]->GetNbinsX();i++){
      double d = hefficiency_d[ifile]->GetBinContent(i+1);
      double n = hefficiency[ifile]->GetBinContent(i+1);
      if (n == 0)
        continue;
      double err = (n/d)*(1-n/d)/d;
      hefficiency[ifile]->SetBinError(i+1,err);
      hefficiency[ifile]->SetBinContent(i+1,n/d);
    }
  }

  hefficiency[0]->SetStats(0);
  hefficiency[0]->SetTitle("");
  hefficiency[0]->GetXaxis()->SetTitle("DOCA to straw center (mm)");
  hefficiency[0]->GetYaxis()->SetTitle("Fractional efficiency");
  hefficiency[0]->SetLineColor(kBlue);
  hefficiency[1]->SetLineColor(kRed);
  hefficiency[0]->SetMarkerStyle(22);
  hefficiency[0]->SetMarkerColor(kBlue);

  lavg[0]->SetLineColor(kBlue);
  lavg[1]->SetLineColor(kRed);
  lavg[0]->SetLineWidth(2);
  lavg[1]->SetLineWidth(2);

  TCanvas *c_eff = new TCanvas("c_eff","c_eff",600,600);
  hefficiency[0]->Draw("E");
//  hefficiency[0]->Draw("hist same");
  hefficiency[1]->Draw("E same");
  hefficiency[1]->Draw("hist same");
  lavg[0]->Draw();
  lavg[1]->Draw();
 
  TLegend *l = new TLegend(0.65,0.65,0.85,0.85);
  l->AddEntry(hefficiency[0],"Data","lepf");
  l->AddEntry(lavg[0],"Data fit","l");
  l->AddEntry(hefficiency[1],"Simulation","lepf");
  l->AddEntry(lavg[1],"Simulation fit","l");
  l->Draw();

  TPaveText *pt = new TPaveText(0.65,0.35,0.85,0.55,"NDC");
  TText *tt = pt->AddText(TString::Format("#epsilon = %.3f #pm %.3f ns",effavg[0],effavgerr[0]));
  tt->SetTextColor(kBlue);
  tt->Draw();
  TText *tt2 = pt->AddText(TString::Format("#epsilon = %.3f #pm %.3f ns",effavg[1],effavgerr[1]));
  tt2->SetTextColor(kRed);
  tt2->Draw();
  pt->Draw();

  theApp.Run();

  return 0;

}

