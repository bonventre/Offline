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
#include "fit/StrawFit.hh"
#include "event/Hit.h"

int main(int argc, char** argv)
{
//  gErrorIgnoreLevel = 1001;

  int strawnum = -1;
  std::string filename_list = "";
  std::string outname = "";
  bool old_data = false;
  std::string help = "./efficiency_data -s <strawnum> -f <filename list> -o <outname> [-u (old data) -p <pmp minimum> -H <fix hdist> -V <fix vdist>]";
  std::string inputcommand = "";
  double fixedhpos = -999;
  double fixedvpos = -999;
  double pmpmin = -9e9;
  for (int i=0;i<argc;i++)
    inputcommand += std::string(argv[i]);
  

  int c;
  while ((c = getopt (argc, argv, "hus:f:o:p:H:V:")) != -1){
    switch (c){
      case 'h':
        std::cout << help << std::endl;
        return 0;
      case 'H':
        fixedhpos = atof(optarg);
        break;
      case 'V':
        fixedvpos = atof(optarg);
        break;
      case 'f':
        filename_list = std::string(optarg);
        break;
      case 's':
        strawnum = atoi(optarg);
        break;
      case 'o':
        outname = std::string(optarg);
        break;
      case 'p':
        pmpmin = atof(optarg);
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

  parse_data_files(filenames, strawnum, alltop_x, allbot_x, alltop_y, allbot_y, has_straw, has_shower, alltimes, alldts, allmaxvals, alltots, old_data, pmpmin);




// find position that gets most of the hits to seed position
  int max_hits = 0;
  double seed_hdist;
  double seed_vdist;
  for (int ih=0;ih<20;ih++){
    for (int iv=0;iv<20;iv++){
      double temp_hdist = ih;
      double temp_vdist = iv;
      int temp_hits = 0;
      for (int i=0;i<alltop_x.size();i++){
        if (has_straw[strawnum][i]){
          if (calculate_DOCA(alltop_x[i],allbot_x[i],temp_hdist,temp_vdist) < 2.5)
            temp_hits++;
//          //FIXME CUT ON MAXVAL TO CUT NOISE
//          if (allmaxvals[strawnum][i] < 80){
//            has_shower[strawnum][i] = true;
//          }
        }
      }
      if (temp_hits > max_hits){
        max_hits = temp_hits;
        seed_hdist = temp_hdist;
        seed_vdist = temp_vdist;
      }
    }
  }
          
  std::cout << "Got " << max_hits << " (" << max_hits/(double) alltop_x.size() << ") at " << seed_hdist << " " << seed_vdist << std::endl;


  std::vector<double> strawbestfit;
  std::vector<double> strawbestfiterrors;
  StrawFit strawfcn(alltop_x,allbot_x,has_straw[strawnum],has_shower[strawnum]);

  std::vector<double> strawseed(4,0);
  strawseed[0] = 0.95;
  strawseed[1] = 0.5;
  strawseed[2] = seed_hdist;
  strawseed[3] = seed_vdist;
  std::vector<double> strawerrors(4,0.1);
  if (fixedhpos > -999){
    strawseed[2] = fixedhpos;
  }
  if (fixedvpos > -999){
    strawseed[3] = fixedvpos;
  }


  ROOT::Minuit2::MnUserParameters strawparams(strawseed,strawerrors);
  ROOT::Minuit2::MnMigrad strawmigrad(strawfcn,strawparams);
  //strawmigrad.Fix((unsigned int) 2);
  //strawmigrad.Fix((unsigned int) 3);
  if (fixedhpos > -999){
    strawmigrad.Fix((unsigned int) 2);
  }
  if (fixedvpos > -999){
    strawmigrad.Fix((unsigned int) 3);
  }

  strawmigrad.SetLimits((unsigned int) 0,0,1);
  strawmigrad.SetLimits((unsigned int) 1,0,3);
  ROOT::Minuit2::FunctionMinimum strawmin = strawmigrad();
  ROOT::Minuit2::MnUserParameters strawresults = strawmin.UserParameters();
  strawbestfit = strawresults.Params();
  strawbestfiterrors = strawresults.Errors();

  double meaneff = strawbestfit[0];
  double exponent = strawbestfit[1];
  double hdist = strawbestfit[2];
  double vdist = strawbestfit[3];


  TFile *fout = new TFile(outname.c_str(),"RECREATE");
  TTree *tout = new TTree("efficiency","efficiency");
  double t_top, t_bot, t_doca, t_pmp;
  int t_hit, t_shower;
  tout->Branch("top",&t_top,"top/D");
  tout->Branch("bot",&t_bot,"bot/D");
  tout->Branch("doca",&t_doca,"doca/D");
  tout->Branch("pmp",&t_pmp,"pmp/D");
  tout->Branch("hit",&t_hit,"hit/I");
  tout->Branch("shower",&t_shower,"shower/I");
  for (int i=0;i<alltop_x.size();i++){
    t_top = alltop_x[i];
    t_bot = allbot_x[i];
    t_doca = calculate_DOCA(alltop_x[i],allbot_x[i],hdist,vdist);
    t_hit = has_straw[strawnum][i];
    t_shower = has_shower[strawnum][i];
    t_pmp = allmaxvals[strawnum][i];
    tout->Fill();
  }
  tout->Write();
  TTree *tresults = new TTree("efficiencyfit","efficiencyfit");
  double t_meaneff, t_exponent, t_hdist, t_vdist;
  double t_meanefferr, t_exponenterr, t_hdisterr, t_vdisterr;
  tresults->Branch("meaneff",&t_meaneff,"meaneff/D");
  tresults->Branch("exponent",&t_exponent,"exponent/D");
  tresults->Branch("hdist",&t_hdist,"hdist/D");
  tresults->Branch("vdist",&t_vdist,"vdist/D");
  tresults->Branch("meanefferr",&t_meanefferr,"meanefferr/D");
  tresults->Branch("exponenterr",&t_exponenterr,"exponenterr/D");
  tresults->Branch("hdisterr",&t_hdisterr,"hdisterr/D");
  tresults->Branch("vdisterr",&t_vdisterr,"vdisterr/D");
  t_meaneff = strawbestfit[0];
  t_exponent = strawbestfit[1];
  t_hdist = strawbestfit[2];
  t_vdist = strawbestfit[3];
  t_meanefferr = strawbestfiterrors[0];
  t_exponenterr = strawbestfiterrors[1];
  t_hdisterr = strawbestfiterrors[2];
  t_vdisterr = strawbestfiterrors[3];

  tresults->Fill();
  tresults->Write();

  TObjString *tl = new TObjString(inputcommand.c_str());
  fout->WriteObject(tl,"log");
  fout->Close();


  std::cout << "Mean efficiency: " << t_meaneff << " +- " << t_meanefferr << std::endl;
  std::cout << "Exponent: " << t_exponent << " +- " << t_exponenterr << std::endl;
  std::cout << "hdist: " << t_hdist << " +- " << t_hdisterr << std::endl;
  std::cout << "vdist: " << t_vdist << " +- " << t_vdisterr << std::endl;

  TH1F *efficiency_straw = new TH1F("eff_straw","eff_straw",100,0,5);
  TH1F *efficiency_straw_d = new TH1F("eff_straw_d","eff_straw",100,0,5);

  TH1F *dt_inside = new TH1F("dt_inside","dt_inside",500,-5,5);
  TH1F *dt_outside = new TH1F("dt_outside","dt_outside",500,-5,5);
  int outside_straw = 0;
  for (int j=0;j<alltop_x.size();j++){
    if (!has_shower[strawnum][j]){
      double straw_dist = calculate_DOCA(alltop_x[j],allbot_x[j],hdist,vdist);
      efficiency_straw_d->Fill(straw_dist);
      if (has_straw[strawnum][j]){
        efficiency_straw->Fill(straw_dist);
        if (straw_dist > 2.5){
          outside_straw++;
   //       std::cout << straw_dist << " : " << allmaxvals[strawnum][j] << " " << alldts[strawnum][j] << std::endl;
          dt_inside->Fill(alldts[strawnum][j]);
        }else{
          dt_outside->Fill(alldts[strawnum][j]);
        }
      }
    }
  }
  std::cout << "Number of hits outside straw: " << outside_straw << std::endl;

  double tot_d = efficiency_straw_d->Integral(1,40);
  double tot_n = efficiency_straw->Integral(1,40);
  std::cout << "AVERAGE: " << tot_n/tot_d << " +- " << 1/tot_d*sqrt(tot_n*(1-tot_n/tot_d)) << std::endl;

  TH1F* hpdf= new TH1F("pdf","pdf",80,0,4);
  for (int j=0;j<hpdf->GetNbinsX();j++){
    double dist = hpdf->GetBinCenter(j+1);
    double pdf_prob = strawbestfit[0] * (1 - 0.5*(1+TMath::Erf((dist-2.5)/(fabs(strawbestfit[1])*sqrt(2)))));
    hpdf->SetBinContent(j+1,pdf_prob);
  }

  //efficiency_straw->Divide(efficiency_straw_d);

  
  for (int i=0;i<efficiency_straw->GetNbinsX();i++){
      double d = efficiency_straw_d->GetBinContent(i+1);
      double n = efficiency_straw->GetBinContent(i+1);
      if (n == 0)
        continue;
      double err = (n/d)*(1-n/d)/d;
      efficiency_straw->SetBinError(i+1,err);
      efficiency_straw->SetBinContent(i+1,n/d);
  }


  efficiency_straw->SetStats(0);
  efficiency_straw->SetTitle("");
  efficiency_straw->GetXaxis()->SetTitle("DOCA to straw center (mm)");
  efficiency_straw->GetYaxis()->SetTitle("Fractional efficiency");
  efficiency_straw->SetLineColor(kBlue);
  efficiency_straw->SetMarkerStyle(22);
  efficiency_straw->SetMarkerColor(kBlue);

  double effavg = tot_n/tot_d;
  double effavgerr = 1/tot_d*sqrt(tot_n*(1-tot_n/tot_d));
  TLine *lavg = new TLine(0,effavg,2.0,effavg);
  lavg->SetLineColor(kBlue);
  lavg->SetLineWidth(2);

  TCanvas *c_eff = new TCanvas("c_eff","c_eff",600,600);
  efficiency_straw->Draw("E");
  hpdf->Draw("hist same");
//  hefficiency[0]->Draw("hist same");
  lavg->Draw();
 
  TLegend *l = new TLegend(0.65,0.65,0.85,0.85);
  l->AddEntry(efficiency_straw,"Data","lepf");
  l->AddEntry(hpdf,"Data fit","l");
  l->AddEntry(lavg,"Average","l");
  l->Draw();

  TPaveText *pt = new TPaveText(0.65,0.35,0.85,0.55,"NDC");
  TText *tt = pt->AddText(TString::Format("#epsilon = %.3f #pm %.3f ns",effavg,effavgerr));
  tt->SetTextColor(kBlue);
  tt->Draw();
  pt->Draw();
  
  //TCanvas *c2 = new TCanvas("c2","c2",600,600);
  //dt_inside->Draw();
  //dt_outside->SetLineColor(kRed);
  //dt_outside->Draw("same");

  theApp.Run();



  return 0;
}
