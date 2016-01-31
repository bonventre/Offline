///////////////////////////////////////////////////////////////////////////////
// class to process calo-digitized waveform
// 
///////////////////////////////////////////////////////////////////////////////
#include "HitMakers/inc/RawWaveformProcess.hh"
#include "art/Framework/Services/Optional/TFileDirectory.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Principal/Run.h"


#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/ConditionsService.hh"
#include "ConditionsService/inc/DAQParams.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"

#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"


namespace mu2e {

//-----------------------------------------------------------------------------
  RawWaveformProcess::RawWaveformProcess(fhicl::ParameterSet const& PSet) :

    WaveformProcess(PSet),
    _debugLevel        (PSet.get<int>   ("debugLevel")),
    _acquisitionEndTime(PSet.get<int>   ("acquisitionEndTime")),
    _digiSampling      (PSet.get<int>   ("digiSampling"))
  {
    _counter        = 0;
    //    _debugHistIndex = 0;
  }


//-----------------------------------------------------------------------------
// destructor
//-----------------------------------------------------------------------------
  RawWaveformProcess::~RawWaveformProcess() {}

  
//------------------------------------------------------------------------------------------//
  void        RawWaveformProcess::processWaveform(double   ADCToMeV, 
						  CaloDigi CaloHit , 
						  RecoCaloDigiCollection &RecoCaloHits) {
    _ADCToMeV  = ADCToMeV;
    
    _amplitude = 0;
    _charge    = 0;

    //reset the fitter
    //    _linearFit->clear();
      
    //lenght of the waveform
    int            wfSize   = CaloHit.waveform().size();

    //time of the first digitized timestamp
    double         t0       = CaloHit.t0();
    double         time, timeBin, content, eDep, chi2;
    
    for (int i=0; i<wfSize; ++i){
      timeBin = i*_digiSampling;
      time    = t0 + timeBin;
      content = CaloHit.waveform().at(i);

      if ( content > _amplitude){
	_amplitude = content;
      }
      _charge += content;
    }

    //convert ADC counts into MeV
    eDep     = _charge*_ADCToMeV*_digiSampling; 
    
    //evaluate the time from the linear fit
    time     =  CaloHit.t0();
    chi2     =  -1;

    //understand how many signals are within the pulse
    _psd  =  _amplitude/_charge;

    RecoCaloHits.push_back( RecoCaloDigi(CaloDigi(CaloHit),
					 eDep,
					 _amplitude, 
					 time, 
					 chi2,
					 _psd    )); 
    

  }
//------------------------------------------------------------------------------------------//

  void        RawWaveformProcess::book(){
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir("CaloDigiDiag");
    
    _hist._hEdep         = tfdir.make<TH1F>("hEdep","Hit energy deposition",200,0.,500);
    _hist._hAmplitude    = tfdir.make<TH1F>("hAmplitude","Waveform amplitude", 2500,0.,2500);

    _hist._hTime         = tfdir.make<TH1F>("hTime","Hit time ", 12000, 0., 2000);

    _hist._hTimeVsE      = tfdir.make<TH2F>("hTimeVsE",
					   "Hit time resolution vs reco charge",
					   500 ,   0., 500,
					   1000, -50.,  50);
    _hist._hTimeVsAmp    = tfdir.make<TH2F>("hTimeVsAmplitude",
					      "Hit amplitude resolution vs reco charge",
					   1250 ,   0., 2500,
					   1000, -20.,  20);

    _hist._hAmpVsE       = tfdir.make<TH2F>("hAmpVsE","Amplitude versus charge; Amplitue [mV]; Energy [MeV]",
					   2500,0.,2500,
					   500 ,   0., 500);
    _hist._hPSDVsAmp     = tfdir.make<TH2F>("hPSDVsAmp","Pile-up disriminator versus amplitude",
					    500 ,   0.,  100,
					    2500,   0., 2500);
    

    _hist._hPSD          = tfdir.make<TH1F>("hPSD","PSD ", 600, 0., 300);
    _hist._hChi2Time     = tfdir.make<TH1F>("hChi2Time","#Chi^2 from time reco alg:#Chi2^2/ndof"   , 100, 0., 10);
    _hist._hNDofTime     = tfdir.make<TH1F>("hNDofTime","nDof from time reconstruction fit; nDof [#]", 200, 0., 10);
    _hist._hDtBinTime    = tfdir.make<TH1F>("hDtBinTime","#Delta t time reco; t_{reco} - t_{bin-edge} [ns]", 400, 0.,_digiSampling);
    
    //histograms for the time reconstruction fit results
    _hist._hFitM       = tfdir.make<TH1F>("hFitM", "Fit - slope; m", 1000, 0, 100);



    _hist._hDt           = tfdir.make<TH1F>("hDt","#Deltat distribution; #Delta t = t_{reco} t_{MC} [ns]", 4000, -100., 100);


    //    int       nDigiSamples =  1965 - _acquisitionEndTime; 
    //    double    mbtime       =  _digiSampling*nDigiSamples;
  
    // for (int i=0; i<20; ++i){
    //   _hist._debugWf[i] = tfdir.make<TH1F>(Form("hDWf%i", i), Form("waveform %i",i),  nDigiSamples, 0, mbtime);
    // }
    //create the  TTree      
    _tree  = tfdir.make<TTree>("Calo", "Calo");
       
    _tree->Branch("counter",      &_counter ,    "counter/I");

    _tree->Branch("time",         &_timeWf  ,        "time/F");
    _tree->Branch("mcTime",       &_mcTime ,      "mcTime/F");
    _tree->Branch("mcMeanTime",   &_mcMeanTime ,  "mcMeanTime/F");
    _tree->Branch("mcEDep",       &_mcEDep ,      "mcEDep/F");

    _tree->Branch("Chi2Time",     &_Chi2Time ,    "Chi2Time/F");
    _tree->Branch("psdWf   ",     &_psdWf      ,  "psdWf/F");

    _tree->Branch("nDof",         &_nDof ,        "nDof/F");
    _tree->Branch("fitM",         &_fitM ,        "fitM/F");
    _tree->Branch("fitQ",         &_fitQ ,        "fitQ/F");

    _tree->Branch("nParticles",   &_nParticles,   "nParticles/I");
    _tree->Branch("pdgId"     ,   &_pdgId     ,   "pdgId[nParticles]/I");

    _tree->Branch("amp",          &_amp ,         "amp/F");
    _tree->Branch("nPeaks",       &_nPeaks ,      "nPeaks/I");
    _tree->Branch("charge",       &_charge ,      "chargexs/F");
  }

  //--------------------------------------------------------------------------------//


  void      RawWaveformProcess::fillDiag(CaloDigiMC*DigiMC, RecoCaloDigiCollection* RecoCaloHits){
    RecoCaloDigi *recoHit;
    int          size = RecoCaloHits->size();

    double       eDep, eDepMax(0);
    
    
    recoHit = &RecoCaloHits->at(size - 1);
    
    _counter     = _hitCounter;
    _psdWf       = _psd;
    _charge      = recoHit->edep();
    
    double     recoTimeBest(0), timeMCBest(0);

    double     recoTime = recoHit->time();
    double     nDof     = 0;
    double     m        = 0;
    double     q        = 0;
    
    _timeWf      = recoTime;
    _Chi2Time    = recoHit->tChi2();
    _nDof        = nDof;
    _amp         = recoHit->amplitude();
    _fitM        = m;
    _fitQ        = q;
    
    double        timeMC     = DigiMC->timeFirst();//meanTime();
    int           nParticles(0);// = DigiMC->nParticles();
    
    _mcTime     = timeMC;
    _mcMeanTime = DigiMC->meanTime();
    _mcEDep     = DigiMC->totalEDep();
    
    
    for (int j=0; j<DigiMC->nParticles(); ++j){
      if (DigiMC->eDep(j) > 1.){
	_pdgId[j] = DigiMC->pdgId(j);
	++nParticles;
      }
      
    }
    _nParticles = nParticles;
    
    _tree->Fill();
    
    ++_hitCounter;


    //now fill the histograms
    _hist._hPSD       ->Fill(_psd);
    _hist._hPSDVsAmp  ->Fill(_psd, recoHit->amplitude());

    _hist._hChi2Time  ->Fill(recoHit->tChi2());
    _hist._hNDofTime  ->Fill(nDof  );
    _hist._hFitM    ->Fill(m   );

    _hist._hEdep      ->Fill(recoHit->edep());
    _hist._hAmplitude ->Fill(recoHit->amplitude());

    _hist._hTime      ->Fill(recoTime);
    _hist._hTimeVsAmp ->Fill(recoTime, recoHit->amplitude());
      

    eDep = recoHit->edep();

    if (eDep > eDepMax){
      eDepMax      = eDep;
      recoTimeBest = recoTime;
      timeMCBest   = timeMC;
    }


      // if (_debugHistIndex < 20){
      // 	double      dt = (recoTime - timeMC);

      // 	if ( ( (_debugHistIndex < 10 ) && (dt < 5.3 && recoHit->edep() > 30.) ) ||
      // 	     ( (_debugHistIndex >= 10) && (dt > 5.7 && recoHit->edep() > 30.) ) ){
      // 	  for (int i=0; i<_wave->GetNbinsX(); ++i){
      // 	    content   = _wave->GetBinContent(i+1);
      // 	    _hist._debugWf[_debugHistIndex] ->SetBinContent(i+1, content);
      // 	    _hist._debugWf[_debugHistIndex] ->SetBinError  (i+1, _wave_point_error);
      // 	  }
      // 	  ++_debugHistIndex;
      // 	}
	
      // }
      //end filling pulses
    
    //    }//end loop on the RecoCaloDigi
    
 
    
    _hist._hDtBinTime ->Fill( (recoTimeBest/_digiSampling - int(recoTimeBest/_digiSampling))*_digiSampling);
    _hist._hDt        ->Fill(recoTimeBest - timeMCBest);
    _hist._hTimeVsE   ->Fill(eDepMax, recoTimeBest- timeMCBest);


  }
  //------------------------------------------------------------------------------------------//

}
