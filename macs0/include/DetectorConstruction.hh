#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"
#include "G4FieldManager.hh"
#include "ElectricFieldSetup.hh"

class G4Box;
class G4Tubs;
class G4UnionSolid;
class G4MultiUnion;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
class DetectorMessenger;
class MagneticField;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:

  DetectorConstruction();
  ~DetectorConstruction();

public:

  virtual G4VPhysicalVolume* Construct();
  void ConstructSDandField();
  G4VPhysicalVolume* macs0(); // MACS version 0, extracted from PRL,82,49

  const
  G4VPhysicalVolume* GetTracker() {return fPhysiTracker;};
  G4double GetTrackerFullLength() {return fTrkLength;};
  G4double GetTargetFullLength()  {return fTgtLength;};
  G4double GetWorldFullLength()   {return fWorldLength;};

  void defineMaterials();
  void SetMagField(G4double);
  void makeSplitTrspTube();
  void makeCombinedTrspTube();
  void makeAccel();

private:
  // -- Materials
  G4Material *fVac, *fAir, *fAl, *fPb, *fSiO2,
    *fXeGas;

  // -- world
  G4Box*             fSolidWorld;
  G4LogicalVolume*   fLogicWorld;
  G4VPhysicalVolume* fPhysiWorld;

  // -- target
  G4Box*             fSolidTarget;
  G4LogicalVolume*   fLogicTarget;
  G4VPhysicalVolume* fPhysiTarget;
  G4Material*        fTargetMater;

  // -- tracker
  G4Tubs*            fSolidTracker;
  G4LogicalVolume*   fLogicTracker;
  G4VPhysicalVolume* fPhysiTracker;

  // -- tracker chambers
  G4Box*             fSolidChamber;
  G4LogicalVolume**  fLogicChamber;
  G4VPhysicalVolume**fPhysiChamber;
  G4Material*        fChamberMater;

  // -- transport tube
  G4MultiUnion*      fSolidTrsp;
  G4LogicalVolume*   fLogicTrsp;
  G4VPhysicalVolume* fPhysiTrsp;
  G4Material*        fTrspMater;
  G4double           fTrspOuterRadius, fTrspLength1, fTrspLength2;
  G4LogicalVolume*   fMagneticLogical;

  // -- accelerator
  G4LogicalVolume*   fLogicAccel;
  G4VPhysicalVolume* fPhysiAccel;


  static G4ThreadLocal MagneticField* fpMagField;
  static G4ThreadLocal G4FieldManager* fpFieldMgr;

  DetectorMessenger* fDetectorMessenger;  // pointer to the Messenger
  G4UserLimits* fStepLimit;            // pointer to user step limits

  G4double fWorldLength;
  G4double fTgtLength;
  G4double fTrkLength, fTrkOuterRadius, fTrkInnerRadius;
  G4int fNbOfChambers;
  G4double fChamberWidth;
  G4double fChamberSpacing;

  static const G4int fMacsTrkNum = 5;
  G4double fMacsTrkLength[fMacsTrkNum], fMacsTrkInnerRadius[fMacsTrkNum],
    fMacsTrkRadialThickness;

  G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

  G4Cache<ElectricFieldSetup*> fEmFieldSetup;

};


#endif
