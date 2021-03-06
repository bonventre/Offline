// Andrei Gaponenko, 2016

#include <string>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Utilities/InputTag.h"

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "CutAndCountAnalysis/inc/CutAndCountAnalysis.hh"

namespace mu2e {
  class CutAndCountFilter : public art::EDFilter {
  public:
    explicit CutAndCountFilter(const fhicl::ParameterSet& pset);
    virtual bool filter(art::Event& event) override;
  private:
    CutAndCountAnalysis an_;
  };

  //================================================================
  CutAndCountFilter::CutAndCountFilter(const fhicl::ParameterSet& pset):
    art::EDFilter{pset},
    an_(pset, *art::ServiceHandle<art::TFileService>())
  {}

  //================================================================
  bool CutAndCountFilter::filter(art::Event& event) {
    return an_.accepted(event);
  }

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CutAndCountFilter);
