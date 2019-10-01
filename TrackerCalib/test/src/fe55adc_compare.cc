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
#include <TF1.h>
#include <TLegend.h>
#include <TApplication.h>

int main(int argc, char** argv)
{
  std::string filename[2];
  filename[0] = "";
  filename[1] = "";
  std::string help = "./compare_fe55adc -d <data fe55 adc fit> -s <sim sddiag>";
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
  
  std::vector<double> pmp[2];
  double t_mean[2], t_meanerr[2], t_width[2], t_widtherr[2];
  TH1F *h[2];
  
  for (int ifile=0;ifile<2;ifile++){
    TFile *f = new TFile(filename[ifile].c_str());
    TTree *t = (TTree*) f->Get("fe55adc");
    double t_pmp;
    t->SetBranchAddress("pmp",&t_pmp);
    for (int i=0;i<t->GetEntries();i++){
      t->GetEntry(i);
      pmp[ifile].push_back(t_pmp);
    }
    TTree *tfit = (TTree*) f->Get("fe55adcfit");
    tfit->SetBranchAddress("peak0mean",&t_mean[ifile]);
    tfit->SetBranchAddress("peak0meanerr",&t_meanerr[ifile]);
    tfit->SetBranchAddress("peak0width",&t_width[ifile]);
    tfit->SetBranchAddress("peak0widtherr",&t_widtherr[ifile]);
    tfit->GetEntry(0);

    h[ifile]= new TH1F(TString::Format("h_%d",ifile),"h",100,0,1000);
    h[ifile]->Sumw2();
    for (int i=0;i<pmp[ifile].size();i++){
      h[ifile]->Fill(pmp[ifile][i]);
    }
    h[ifile]->Scale(1.0/h[ifile]->Integral());
    h[ifile]->SetDirectory(0);
  }

  std::cout << "Sim Peak: " << t_mean[1] << " +- " << t_meanerr[1] << std::endl;
  std::cout << "Sim Width: " << t_width[1] << " +- " << t_widtherr[1] << std::endl;
  std::cout << "Data Peak: " << t_mean[0] << " +- " << t_meanerr[0] << std::endl;
  std::cout << "Data Width: " << t_width[0] << " +- " << t_widtherr[0] << std::endl;
 
  h[0]->SetStats(0);
  h[0]->SetTitle("");
  h[0]->GetXaxis()->SetTitle("ADC Peak minus pedestal (ADC counts)");
  h[0]->GetYaxis()->SetTitle("");
  h[0]->SetLineColor(kBlue);
  h[1]->SetLineColor(kRed);
  h[0]->SetMarkerStyle(22);
  h[0]->SetMarkerColor(kBlue);

  h[0]->Draw("E");
  h[1]->Draw("E same");
  h[1]->Draw("HIST same");

  TLegend *l = new TLegend(0.65,0.65,0.85,0.85);
  l->AddEntry(h[0],"Data","lepf");
  l->AddEntry(h[1],"Simulation","lepf");
  l->Draw();

  theApp.Run();

  return 0;
}
