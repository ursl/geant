#ifndef B1DetectorConstruction_h
#define B1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4AssemblyVolume; 

class B1DetectorConstruction : public G4VUserDetectorConstruction {
public:
  B1DetectorConstruction();
  virtual ~B1DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();

  void placeFEE();
  G4AssemblyVolume*  makeFEE(G4LogicalVolume *);
  
  G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }
  
protected:
  G4LogicalVolume*  fVolume;
  G4LogicalVolume*  fScoringVolume;
};

#endif

