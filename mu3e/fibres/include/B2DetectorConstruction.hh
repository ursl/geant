#ifndef B2DetectorConstruction_h
#define B2DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4AssemblyVolume; 

class B2DetectorConstruction : public G4VUserDetectorConstruction {
public:
  B2DetectorConstruction();
  virtual ~B2DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();

  void placeSMB();
  G4AssemblyVolume*  makeSMB(G4LogicalVolume *);
  
  G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
  
protected:
  G4LogicalVolume*  fVolume;
  G4LogicalVolume*  fScoringVolume;
};

#endif

