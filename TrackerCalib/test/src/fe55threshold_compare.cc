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
#include <TH1F.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TApplication.h>

int main(int argc, char** argv)
{
  std::string filename[2];
  filename[0] = "";
  filename[1] = "";
  std::string help = "./fe55threshold_compare -d <data fe55 threshold fit> -s <sim fe55 threshold fit>";
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

  if (filename[0].size() == 0 || filename[1].size() == 0){
    std::cout << help << std::endl;
    return 1;
  }

  int argc2 = 0;char **argv2;TApplication theApp("tapp", &argc2, argv2);
  
  std::vector<double> measured_thresh[2];
  std::vector<double> thresh[2];
  std::vector<double> width[2];
  double t_threshoffset, t_threshslope, t_threshoffseterr, t_threshslopeerr;
  double t_widthoffset, t_widthslope, t_widthoffseterr, t_widthslopeerr;
  TGraph *gthresh[2];
  TGraph *gwidth[2];
  TF1 *fthresh = new TF1("fthresh","[0] + [1]*x",0,50);
  TF1 *fwidth = new TF1("fwidth","[0] + [1]*x",0,50);
  
  for (int ifile=0;ifile<2;ifile++){
    TFile *f = new TFile(filename[ifile].c_str());
    TTree *t = (TTree*) f->Get("fe55threshold");
    double t_measured_thresh;
    double t_thresh;
    double t_width;
    t->SetBranchAddress("measured_thresh",&t_measured_thresh);
    t->SetBranchAddress("thresh",&t_thresh);
    t->SetBranchAddress("width",&t_width);
    for (int i=0;i<t->GetEntries();i++){
      t->GetEntry(i);
      measured_thresh[ifile].push_back(t_measured_thresh);
      thresh[ifile].push_back(t_thresh);
      width[ifile].push_back(t_width);
    }
    if (ifile == 0){
      TTree *tfit = (TTree*) f->Get("fe55thresholdfit");
      tfit->SetBranchAddress("threshoffset",&t_threshoffset);
      tfit->SetBranchAddress("threshslope",&t_threshslope);
      tfit->SetBranchAddress("threshoffseterr",&t_threshoffseterr);
      tfit->SetBranchAddress("threshslopeerr",&t_threshslopeerr);
      tfit->SetBranchAddress("widthoffset",&t_widthoffset);
      tfit->SetBranchAddress("widthslope",&t_widthslope);
      tfit->SetBranchAddress("widthoffseterr",&t_widthoffseterr);
      tfit->SetBranchAddress("widthslopeerr",&t_widthslopeerr);
      tfit->GetEntry(0);
      fthresh->SetParameter((unsigned)0,t_threshoffset);
      fthresh->SetParameter((unsigned)1,t_threshslope);
      fwidth->SetParameter((unsigned)0,t_widthoffset);
      fwidth->SetParameter((unsigned)1,t_widthslope);
    }

    gthresh[ifile] = new TGraph(measured_thresh[ifile].size(),&measured_thresh[ifile][0],&thresh[ifile][0]);
    gwidth[ifile] = new TGraph(measured_thresh[ifile].size(),&measured_thresh[ifile][0],&width[ifile][0]);

//    gthresh[ifile]->SetDirectory(0);
//    gwidth[ifile]->SetDirectory(0);
  }

  gthresh[0]->SetMarkerStyle(21);
  gthresh[0]->SetMarkerColor(kBlack);
  gthresh[0]->SetLineColor(kBlack);
  gthresh[1]->SetMarkerStyle(22);
  gthresh[1]->SetMarkerColor(kBlue);
  gthresh[1]->SetLineColor(kBlue);

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TMultiGraph *mg = new TMultiGraph("mg","mg");
  mg->Add(gthresh[0]);
  mg->Add(gthresh[1]);
  mg->SetTitle("");
  mg->Draw("ALP");
  mg->GetXaxis()->SetTitle("Set threshold (mV)");
  mg->GetYaxis()->SetTitle("Output signal threshold (peak minus pedestal counts)");
  //gthresh[1]->Draw("same LP");
  fthresh->SetLineColor(kRed);
  fthresh->Draw("same");

  TLegend *l = new TLegend(0.65,0.65,0.85,0.85);
  l->AddEntry(gthresh[0],"Data","lp");
  l->AddEntry(gthresh[1],"Simulation","lp");
  l->AddEntry(fthresh,"Data fit","l");
  l->Draw();

  gwidth[0]->SetMarkerStyle(21);
  gwidth[0]->SetMarkerColor(kBlack);
  gwidth[0]->SetLineColor(kBlack);
  gwidth[1]->SetMarkerStyle(22);
  gwidth[1]->SetMarkerColor(kBlue);
  gwidth[1]->SetLineColor(kBlue);

  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  TMultiGraph *mg2 = new TMultiGraph("mg2","mg2");
  mg2->Add(gwidth[0]);
  mg2->Add(gwidth[1]);
  mg2->SetTitle("");
  mg2->Draw("ALP");
  mg2->GetXaxis()->SetTitle("Set threshold (mV)");
  mg2->GetYaxis()->SetTitle("Signal cutoff width (ADC counts)");
  //gthresh[1]->Draw("same LP");
  fwidth->SetLineColor(kRed);
  fwidth->Draw("same");

  TLegend *l2 = new TLegend(0.65,0.65,0.85,0.85);
  l2->AddEntry(gwidth[0],"Data","lp");
  l2->AddEntry(gwidth[1],"Simulation","lp");
  l2->AddEntry(fwidth,"Data fit","l");
  l2->Draw();



  theApp.Run();

  return 0;
}
