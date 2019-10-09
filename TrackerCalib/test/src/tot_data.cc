#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>

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
#include <TLegend.h>
#include <TPaveText.h>
#include <TText.h>
#include <TGraphAsymmErrors.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TCanvas.h>

#include "prototype/Utils.hh"
#include "prototype/Parser.hh"
#include "event/Hit.h"

int main(int argc, char** argv)
{
//  gErrorIgnoreLevel = 1001;

  int strawnum = -1;
  std::string filename_list = "";
  std::string outname = "";
  double timeoffset = 0;
  bool old_data = false;
  std::string help = "./tot_data -s <strawnum> -f <filename list> -o <outname> [-u (old data)]";
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
          std::cout << "Unknown option `-" << (char) optopt << "'." << std::endl;
        return 1;
    }
  }

  if (strawnum == -1 || filename_list.size() == 0 || outname.size() == 0){
    std::cout << help << std::endl;
    return 1;
  }
  
  int argc2 = 0;char **argv2;TApplication theApp("tapp", &argc2, argv2);

  std::vector<std::string> filenames;
  std::string line;
  std::ifstream fnlist(filename_list.c_str());
  while (std::getline(fnlist,line)){
    filenames.push_back(line);
  }

  std::vector<double> alltop_x; // all hits top pixel avg perpendicular to straw [pixel number (50um)]
  std::vector<double> allbot_x; // all hits bot pixel avg perpendicular to straw [pixel number (50um)]
  std::vector<double> alltop_y; // all hits top pixel avg along straw [pixel number (250um)]
  std::vector<double> allbot_y; // all hits bot pixel avg along straw [pixel number (250um)]

  std::vector<std::vector<bool> > has_straw(8, std::vector<bool>()); // has a hit in straw i
  std::vector<std::vector<bool> > has_shower(8, std::vector<bool>()); // has a hit in a straw not neighbors with straw i

  std::vector<std::vector<double> > alltimes(8, std::vector<double>()); // drift time measurement in straw i
  std::vector<std::vector<double> > alldts(8, std::vector<double>()); // NOT USED
  std::vector<std::vector<double> > allmaxvals(8, std::vector<double>()); // NOT USED
  std::vector<std::vector<double> > alltots(8, std::vector<double>()); // NOT USED

  parse_data_files(filenames, strawnum, alltop_x, allbot_x, alltop_y, allbot_y, has_straw, has_shower, alltimes, alldts, allmaxvals, alltots, old_data);

  double hdist, vdist, pixelseparation;

  TFile *f = new TFile("dataResults/new_1425_12mv_10mm__transverse_straw3.root");
  TTree *tfit = (TTree*) f->Get("transversefit");
  tfit->SetBranchAddress("hdist",&hdist);
  tfit->SetBranchAddress("vdist",&vdist);
  tfit->SetBranchAddress("pixelseparation",&pixelseparation);
  tfit->SetBranchAddress("timeoffset",&timeoffset);
  tfit->GetEntry(0);


  //TH2F *h1 = new TH2F("h1","h1",32,0,64,120,-20,100);
  TH2F *h1 = new TH2F("h1","h1",32,0,64,100,0,5);
  TH1F *h1d = new TH1F("h1d","h1d",100,-5,5);
  TH1F *h2d = new TH1F("h1d","h1d",100,-5,5);
//  timeoffset = -20;

  for (int i=0;i<alltimes[strawnum].size();i++){
    if (has_straw[strawnum][i] && !has_shower[strawnum][i]){
      double doca = calculate_DOCA(alltop_x[i],allbot_x[i],hdist,vdist,pixelseparation);
      double x,y;
      calculate_relative_position(alltop_x[i], allbot_x[i], hdist, vdist, pixelseparation, x, y);
      //if (alltimes[strawnum][i]-timeoffset < -20 || alltimes[strawnum][i]-timeoffset > 20)
      //  continue;
      //if (alltots[strawnum][i]+alltimes[strawnum][i]-timeoffset > 45)
      //  h1->Fill(alltots[strawnum][i],alltimes[strawnum][i]-timeoffset);
      //if (alltots[strawnum][i]+doca*16 > 42.5)
      h1->Fill(alltots[strawnum][i],doca);
      //h1->Fill(y,alltots[strawnum][i]+alltimes[strawnum][i]-timeoffset);
      //h1d->Fill(alltots[strawnum][i]+alltimes[strawnum][i]-timeoffset);
      h1d->Fill(doca-(40-alltots[strawnum][i])/16.);
      h2d->Fill(doca-(50-alltots[strawnum][i])/16.);
    }
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
  hp->Draw();
//  h1->Draw();
//  h1d->Draw();
//  h2d->SetLineColor(kRed);
//  h2d->Draw("same");

  theApp.Run();



  return 0;
}
