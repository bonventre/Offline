//Author: S Middleton
//Purpose: Fit cosmic tracks within the tracker

#ifndef CosmicReco_CosmicTrackFit_HH
#define CosmicReco_CosmicTrackFit_HH

#ifndef __GCCXML__
#include "fhiclcpp/ParameterSet.h"
#endif/*__GCCXML__*/

//Mu2e Cosmics:
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"

// Products
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"

//Drift:
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerConditions/inc/StrawPhysics.hh"
#include "TrackerConditions/inc/StrawDrift.hh"

// Math
#include "Math/VectorUtil.h"
#include "Math/Vector2D.h"

//C++
#include <vector>
#include <utility>
#include <string>
#include <cmath>
#include <algorithm>

// Framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

//ROOT
#include "TMatrixD.h"

namespace mu2e 
{
    class Tracker;
    class CosmicTrackFit
    {
     public:	

	struct Config{
	      using Name=fhicl::Name;
	      using Comment=fhicl::Comment;
	      fhicl::Atom<string> dontuseflag {Name("DoNotUseFlag"),Comment("if set to OK then save the track"), "Outlier"};
              fhicl::Atom<unsigned> minnsh {Name("minNStrawHits"), Comment("minimum number of straw hits "),1};
	      fhicl::Atom<unsigned> n_outliers{Name("Noutliers"),Comment("maximum number of outliers allowed in track fit"),2};
    	      fhicl::Atom<unsigned> maxniter{Name("maxNiter"), Comment("Maximum allowed number of iterations before considered unconverged in seed fit"),10};
    	      
	      fhicl::Atom<float> maxd{Name("MaxTrackLength"),Comment("The maxiumum allowed length of track") ,2000.};
              fhicl::Atom<bool> UseTSeedDirection{Name("UseTSeedDirection"),Comment("Uses the direction in the input tseed to initialize fit"), false};
    	};
		
		explicit CosmicTrackFit(const Config& conf);
    		virtual ~CosmicTrackFit(){};

                XYZVec ConvertPointToDetFrame(XYZVec vec);

                void Chi2Fit(const char* title, CosmicTrackSeed &tseed, art::Event const& event, ComboHitCollection const& chcol, std::vector<StrawHitIndex> &panelHitIdxs);
		TrackEquation ConvertFitToDetectorFrame(TrackAxes axes, XYZVec Position, XYZVec Direction, bool isseed, bool det);
		
                bool goodTrack(CosmicTrack& track);
		
                const Tracker*            _tracker;
    		void  setTracker    (const Tracker*    Tracker) { _tracker     = Tracker; }
                
		bool use_hit(ComboHit const&) const;
  		bool use_track(double length) const;
	private:
		Config _conf;
  		
                StrawHitFlag _dontuseflag; //some flags for removn
    		unsigned _minnsh;  // minimum # of StrawHits
		unsigned _n_outliers; //number of significant outliers/number of hits in track..
    		unsigned _maxniter; // maxium # of iterations to global minimum         
		float _maxd;//unused
                bool _useTSeedDirection;
  };
	

}
#endif
