#ifndef VOLUMEINFO_HH
#define VOLUMEINFO_HH
//
// Information about a physical volume.  Used by Mu2eWorld and its utility routines.
// The center information is not fully general: it does not know about rotations
// and is useful only for the top few levels of the detector.
// 
// 
// $Id: VolumeInfo.hh,v 1.4 2010/08/19 15:46:07 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/08/19 15:46:07 $
//
// Original author Rob Kutschke
//

#include <string>

#include "G4ThreeVector.hh"

class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;

namespace mu2e {

  class VolumeInfo{

  public:

    // Depracated: this will go away.
    VolumeInfo():
      name(),
      solid(0),
      logical(0),
      physical(0),
      centerInParent(),
      centerInWorld(){}

    VolumeInfo( const std::string&       pName,
                const CLHEP::Hep3Vector& inParent,
                const CLHEP::Hep3Vector& parentInWorld);

    // Compiler written versions will be correct for:
    // destructor, copy constructor, assignment operator.

    // The name of this volume as known to G4.
    std::string name;

    // Non-owning pointers to volume information.
    G4VSolid*          solid;
    G4LogicalVolume*   logical;
    G4VPhysicalVolume* physical;
    
    // Location information in two coordinate systems.
    G4ThreeVector      centerInParent;
    G4ThreeVector      centerInWorld;

    G4ThreeVector centerInMu2e() const { return centerInWorld-_Mu2eOriginInWorld; }

    static void setMu2eOriginInWorld( const G4ThreeVector& origin ){
      _Mu2eOriginInWorld = origin;
    }

  private:

    static G4ThreeVector _Mu2eOriginInWorld;

  };

}

#endif
