#ifndef RecoDataProducts_CrvRecoPulses_hh
#define RecoDataProducts_CrvRecoPulses_hh
//
// $Id: $
// $Author: ehrlich $
// $Date: 2014/08/07 01:33:41 $
//
// Contact person Ralf Ehrlich
//

#include <vector>

namespace mu2e 
{
  class CrvRecoPulses
  {
    public:

    CrvRecoPulses() {}

    struct CrvSingleRecoPulse
    {
      int    _PEs;
      double _leadingEdge;
      double _pulseHeight;
      CrvSingleRecoPulse(int PEs, double leadingEdge, double pulseHeight) : _PEs(PEs), 
                                                                            _leadingEdge(leadingEdge), 
                                                                            _pulseHeight(pulseHeight) {}
      CrvSingleRecoPulse();
    };

    std::vector<CrvSingleRecoPulse> &GetRecoPulses(int fiberNumber, int side);
    std::vector<CrvSingleRecoPulse> &GetRecoPulses(int SiPMNumber); 

    const std::vector<CrvSingleRecoPulse> &GetRecoPulses(int fiberNumber, int side) const;
    const std::vector<CrvSingleRecoPulse> &GetRecoPulses(int SiPMNumber) const;

    unsigned int GetNumberOfRecoPulses(int fiberNumber, int side);
    unsigned int GetNumberOfRecoPulses(int SiPMNumber);

    private:

    static int  FindSiPMNumber(int fiberNumber, int side);
    static void CheckSiPMNumber(int SiPMNumber);

    std::vector<CrvSingleRecoPulse> _crvPulses[4];
  };
}

#endif /* RecoDataProducts_CrvRecoPulses_hh */