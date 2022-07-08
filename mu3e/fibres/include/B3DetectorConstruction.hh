#ifndef B3DetectorConstruction_h
#define B3DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Polycone.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4AssemblyVolume; 

class B3DetectorConstruction : public G4VUserDetectorConstruction {
public:
  B3DetectorConstruction();
  virtual ~B3DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();

  G4Polycone*  makeColdFlange(G4LogicalVolume *);
  
  G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
  
protected:
  G4LogicalVolume*  fVolume;
  G4LogicalVolume*  fScoringVolume;
};

#endif

