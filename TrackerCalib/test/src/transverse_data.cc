#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <numeric>

#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>

#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>


#include "prototype/Utils.hh"
#include "prototype/Parser.hh"
#include "fit/DataFit.hh"
#include "fit/TrackerCalibFit.hh"
#include "event/Hit.h"

#define SEED_SIGMA 1
#define SEED_TAU 10

bool compare(const std::pair<double, double>&i, const std::pair<double, double>&j)
{
      return i.first < j.first;
}

int main(int argc, char** argv)
{
  float voltage = 1400;
  int strawnum = -1;
  std::string filename_list = "";
  std::string outname = "";
  bool old_data = false;
  std::string help = "./transverse_data -s <strawnum> -f <filename list> -o <outname> [-v voltage] [-u (old data)]";
  std::string inputcommand = "";
  for (int i=0;i<argc;i++)
    inputcommand += std::string(argv[i]);
  

  int c;
  while ((c = getopt (argc, argv, "hus:f:o:v:")) != -1){
    switch (c){
      case 'h':
        std::cout << help << std::endl;
        return 0;
      case 'v':
        voltage = atof(optarg);      
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
      case 'u':
        old_data = true;
        break;
      case '?':
        if (optopt == 'v' || optopt == 'f' || optopt == 's' || optopt == 'o')
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

  std::vector<double> tops;
  std::vector<double> bots;
  std::vector<double> times;
  std::vector<double> pmps;
  
  for (int j=0;j<alltop_x.size();j++){
    if (has_straw[strawnum][j] && !has_shower[strawnum][j]){
      tops.push_back(alltop_x[j]);
      bots.push_back(allbot_x[j]);
      times.push_back(alltimes[strawnum][j]);
      pmps.push_back(allmaxvals[strawnum][j]);
    }
  }
  std::cout << times.size() << " events for straw " << strawnum << std::endl;


  // find timing rising edge to seed time offset
  double seed_time = 0;
  TH1F *h = new TH1F("","",2000,-100,1900);
  for (int i=0;i<times.size();i++){
    h->Fill(times[i]);
  }
  for (int i=0;i<h->GetNbinsX();i++){
    if (h->GetBinContent(i+1) > h->GetMaximum()/2.){
      std::cout << "SEED AT " << h->GetBinCenter(i+1) << std::endl;
      seed_time = h->GetBinCenter(i+1);
      break;
    }
  }

  // find position that gets most of the hits to seed position
  int max_hits = 0;
  double seed_hdist;
  double seed_vdist;
  for (int ih=0;ih<20;ih++){
    for (int iv=0;iv<20;iv++){
      double temp_hdist = ih;
      double temp_vdist = iv;
      int temp_hits = 0;
      for (int i=0;i<times.size();i++){
        if (calculate_DOCA(tops[i],bots[i],temp_hdist,temp_vdist) < 2.5)
          temp_hits++;
      }
      if (temp_hits > max_hits){
        max_hits = temp_hits;
        seed_hdist = temp_hdist;
        seed_vdist = temp_vdist;
      }
    }
  }
          
  std::cout << "Got " << max_hits << " (" << max_hits/(double) times.size() << ") at " << seed_hdist << " " << seed_vdist << std::endl;


  // Set up seeds for fit
  std::vector<double> seeds(7,0);
  std::vector<double> errors(7,0);
  seeds[0] = SEED_TAU;
  seeds[1] = SEED_SIGMA;
  seeds[2] = seed_time;
  seeds[3] = 0;
  seeds[4] = seed_hdist;
  seeds[5] = seed_vdist;
  seeds[6] = 20.2;

  errors[0] = 0.1;
  errors[1] = 0.1;
  errors[2] = 1;
  errors[3] = 0;
  errors[4] = 0.1;
  errors[5] = 0.1;
  errors[6] = 0.1;

  std::vector<double> constraint_means(7,0);
  std::vector<double> constraints(7,0);
  constraint_means[6] = 20.2;
  constraints[6] = 0.2;

  // fit using MCMC
  std::vector<std::vector<double> > mcmcresults;
  std::vector<double> mcmcnlls;
  DataFit fit(tops,bots,times,constraint_means,constraints,1,voltage);

  int nsteps = 50000;
  fit.metropolis(seeds,errors,mcmcresults,mcmcnlls,nsteps);

  // extract the best fit value from the results
  std::vector<double> mcmcbestfit;
  double mcmcnll = 1e10;
  std::vector<std::vector<double> > jump_vectors;
  std::vector<std::vector<std::pair<double, double> > > valnll_pair;
  std::vector<double> stddev_errors(7,0);
  for (int j=0;j<7;j++){
    jump_vectors.push_back(std::vector<double>());
    valnll_pair.push_back(std::vector<std::pair<double, double> >());
  }
  for (int i=0;i<mcmcresults.size();i++){
    if (mcmcnlls[i] < mcmcnll){
      mcmcnll = mcmcnlls[i];
      mcmcbestfit = mcmcresults[i];
    }
    for (int j=0;j<7;j++){
      jump_vectors[j].push_back(mcmcresults[i][j]);
      valnll_pair[j].push_back(std::make_pair(mcmcresults[i][j],mcmcnlls[i]));
    }
  }


  // calculate the standard deviation for each parameter
  // and create a histogram of the 1-D projections
  std::vector<TH1F*> jump_hists;
  std::vector<std::vector<TH2F*> > jump_hists2d;
  for (int j=0;j<7;j++){
    double sum = std::accumulate(jump_vectors[j].begin(), jump_vectors[j].end(), 0.0);
    double mean = sum / jump_vectors[j].size();
    double sq_sum = std::inner_product(jump_vectors[j].begin(), jump_vectors[j].end(), jump_vectors[j].begin(), 0.0);
    stddev_errors[j] = std::sqrt(sq_sum / jump_vectors[j].size() - mean * mean);
    if (stddev_errors[j] <= 0)
      stddev_errors[j] = 1;

    jump_hists.push_back(new TH1F("","",100,mcmcbestfit[j]-stddev_errors[j]*3,mcmcbestfit[j]+stddev_errors[j]*3));
    jump_hists2d.push_back(std::vector<TH2F*>());
  }
  for (int j=0;j<7;j++){
    for (int k=0;k<7;k++){
      jump_hists2d[j].push_back(new TH2F("","",100,mcmcbestfit[j]-stddev_errors[j]*3,mcmcbestfit[j]+stddev_errors[j]*3,100,mcmcbestfit[k]-stddev_errors[k]*3,mcmcbestfit[k]+stddev_errors[k]*3));
    }
  }
      
  for (int i=0;i<mcmcresults.size();i++){
    for (int j=0;j<7;j++){
      jump_hists[j]->Fill(mcmcresults[i][j]);
      for (int k=0;k<7;k++){
        jump_hists2d[j][k]->Fill(mcmcresults[i][j],mcmcresults[i][k]);
      }
    }
  }


  std::cout << "MCMC NLL: " << mcmcnll << std::endl;
  // calculate uncertainties as 66% of posterior and as delta nll of 0.5
  double minbayesian[7];
  double maxbayesian[7];
  double minfrequentist[7];
  double maxfrequentist[7];
  for (int j=0;j<7;j++){
    std::sort(valnll_pair[j].begin(),valnll_pair[j].end(),compare);
    minbayesian[j] = 9e9;
    maxbayesian[j] = -9e9;
    minfrequentist[j] = 9e9;
    maxfrequentist[j] = -9e9;
    int total_steps = 0;
    for (int i=valnll_pair[j].size()/2;i<valnll_pair[j].size();i++){
      if (total_steps < 0.6827*nsteps/2)
        maxbayesian[j] = valnll_pair[j][i].first;
      if (valnll_pair[j][i].second < mcmcnll + 0.5)
        maxfrequentist[j] = valnll_pair[j][i].first;
      total_steps++;
    }
    total_steps = 0;
    for (int i=valnll_pair[j].size()/2;i>=0;i--){
      if (total_steps < 0.6827*nsteps/2)
        minbayesian[j] = valnll_pair[j][i].first;
      if (valnll_pair[j][i].second < mcmcnll + 0.5)
        minfrequentist[j] = valnll_pair[j][i].first;
      total_steps++;
    }
  }
  std::vector<std::string> names;
  names.push_back("tau (ns)");
  names.push_back("sigma (ns)");
  names.push_back("time offset (ns)");
  names.push_back("background (frac)");
  names.push_back("hdist (mm)");
  names.push_back("vdist (mm)");
  names.push_back("pixel separation (mm)");


  for (int i=0;i<7;i++){

    std::cout << i << names[i] << " : " << mcmcbestfit[i] << " +- " << stddev_errors[i] << " (" << jump_hists[i]->GetRMS() << ") +" << (maxbayesian[i]-mcmcbestfit[i]) << " -" << (mcmcbestfit[i]-minbayesian[i]) << " +" << (maxfrequentist[i]-mcmcbestfit[i]) << " -" << (mcmcbestfit[i]-minfrequentist[i]) << std::endl;
  }
//  int argc2 = 0;char **argv2;TApplication theApp("tapp", &argc2, argv2);

  
  TFile *fout = new TFile(outname.c_str(),"RECREATE");
  TTree *tout = new TTree("transverse","transverse");
  double t_time;
  double t_doca;
  double t_pmp;
  tout->Branch("time",&t_time,"time/D");
  tout->Branch("doca",&t_doca,"doca/D");
  tout->Branch("pmp",&t_pmp,"pmp/D");
  for (int i=0;i<times.size();i++){
    double doca = calculate_DOCA(tops[i],bots[i],mcmcbestfit[4],mcmcbestfit[5],mcmcbestfit[6]);
    t_time = times[i];
    t_doca = doca;
    t_pmp = pmps[i];
    tout->Fill();
  }
  tout->Write();
  TTree *tresults = new TTree("transversefit","transversefit");
  double t_tau, t_sigma, t_timeoffset, t_background, t_hdist, t_vdist, t_pixelseparation;
  double t_tauerr, t_sigmaerr, t_timeoffseterr, t_backgrounderr, t_hdisterr, t_vdisterr, t_pixelseparationerr;
  double t_nll, t_voltage;
  tresults->Branch("nll",&t_nll,"nll/D");
  tresults->Branch("voltage",&t_voltage,"voltage/D");
  tresults->Branch("tau",&t_tau,"tau/D");
  tresults->Branch("sigma",&t_sigma,"sigma/D");
  tresults->Branch("timeoffset",&t_timeoffset,"timeoffset/D");
  tresults->Branch("background",&t_background,"background/D");
  tresults->Branch("hdist",&t_hdist,"hdist/D");
  tresults->Branch("vdist",&t_vdist,"vdist/D");
  tresults->Branch("pixelseparation",&t_pixelseparation,"pixelseparation/D");
  tresults->Branch("tauerr",&t_tauerr,"tauerr/D");
  tresults->Branch("sigmaerr",&t_sigmaerr,"sigmaerr/D");
  tresults->Branch("timeoffseterr",&t_timeoffseterr,"timeoffseterr/D");
  tresults->Branch("backgrounderr",&t_backgrounderr,"backgrounderr/D");
  tresults->Branch("hdisterr",&t_hdisterr,"hdisterr/D");
  tresults->Branch("vdisterr",&t_vdisterr,"vdisterr/D");
  tresults->Branch("pixelseparationerr",&t_pixelseparationerr,"pixelseparationerr/D");
  t_nll = mcmcnll;
  t_voltage = voltage;
  t_tau = mcmcbestfit[0];
  t_sigma = mcmcbestfit[1];
  t_timeoffset = mcmcbestfit[2];
  t_background = mcmcbestfit[3];
  t_hdist = mcmcbestfit[4];
  t_vdist = mcmcbestfit[5];
  t_pixelseparation = mcmcbestfit[6];
  t_tauerr = stddev_errors[0];
  t_sigmaerr = stddev_errors[1];
  t_timeoffseterr = stddev_errors[2];
  t_backgrounderr = stddev_errors[3];
  t_hdisterr = stddev_errors[4];
  t_vdisterr = stddev_errors[5];
  t_pixelseparationerr = stddev_errors[6];
  tresults->Fill();
  tresults->Write();

  TObjString *tl = new TObjString(inputcommand.c_str());
  fout->WriteObject(tl,"log");

  TCanvas *cslice = new TCanvas("cslice","cslice",800,900);
  cslice->Divide(2,3);
  for (int j=0;j<7;j++){
    if (j < 3)
      cslice->cd(j+1);
    if (j == 3)
      continue;
    if (j > 3)
      cslice->cd(j);
    jump_hists[j]->Draw();
  }
  cslice->Write();
  for (int k=0;k<7;k++){
    if (k == 3)
      continue;
    TCanvas *cslicek = new TCanvas(TString::Format("cslice_%d",k),"cslice",800,900);
    cslicek->Divide(2,3);
    for (int j=0;j<7;j++){
      if (j < 3)
        cslicek->cd(j+1);
      if (j == 3)
        continue;
      if (j > 3)
        cslicek->cd(j);
      jump_hists2d[k][j]->Draw();
    }
    cslicek->Write();
  }
  fout->Close();


  return 0;
}

