#include "Mu2eUtilities/inc/SimParticleGetTau.hh"

// Framework includes
#include "art/Persistency/Common/Ptr.h"

// cetlib includes
#include "cetlib/exception.h"

namespace mu2e {
  
  //==========================================================================
  double SimParticleGetTau::calculate( const art::Ptr<SimParticle>& p,
                                       const VspMC& hitColls,
                                       const std::vector<int>& decayOffCodes,
                                       const PhysicsParams& gc ){
    
    double tau = p->endProperTime() / gc.getParticleLifetime(p->pdgId());
    
    // The mu2ePrimary code means that G4 track was created by our PrimaryGeneratorAction.
    // If the particle is "mu2ePrimary" but still has a parent, it is a continuation
    // of a particle from the previous simulation stage, and their proper times should
    // be combined.

    art::Ptr<SimParticle> part (p);
    while(part->parent().isNonnull()) {

      if((part->creationCode() == ProcessCode::mu2ePrimary)) {

        // The current particle is a continuation from the previous stage,
        // not a physically different particle.  We need to find its record
        // in the StepPointMC collections to get the correct proper time.
        part = part->parent();

        // Find matches in hit collections
        unsigned counter (0);
        const StepPointMC* spMC(nullptr);
        for ( const auto& hitColl : hitColls ) {
          std::for_each( hitColl.begin(), hitColl.end(),
                         [&](const StepPointMC& sp){
                           if ( sp.simParticle().key() == part.key() ) {
                             spMC = &sp;
                             counter++;
                           }
                         } );

        }

        if      ( counter == 0 ) throw cet::exception("StepPointMC") << " Non-existent StepPointMC-SimParticle assignment! " ;
        else if ( counter  > 1 ) throw cet::exception("StepPointMC") << " Ambiguous StepPointMC-SimParticle assignment! " ;
        else  tau += spMC->properTime() / gc.getParticleLifetime(part->pdgId());
      }
      else {
        // The current particle was produced by a G4 physics process.
        // See if proper time of its ancestor should be included.
        part = part->parent();
        if ( std::binary_search( decayOffCodes.begin(), decayOffCodes.end(), int(part->pdgId()) ) ) {
          tau += part->endProperTime() / gc.getParticleLifetime(part->pdgId());
        }
      }
    } // loop up to the primary
    
    return tau;
  }
  
} // namespace mu2e