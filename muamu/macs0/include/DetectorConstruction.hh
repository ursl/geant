#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"
#include "G4FieldManager.hh"
//#include "ElectricFieldSetup.hh"

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
//class ElectricFieldSetup;
class F02ElectricFieldSetup;

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
  G4double GetTargetFullLength()  {return fTargetLength;};
  G4double GetWorldFullLength()   {return fWorldLength;};

  void defineMaterials();
  void parseCmd(G4String cmd);
  void setVerbose(int v) {fVerbose = v;}

  // -- communications methods
  void SetMagField(G4double);
  void SetTargetMaterial(G4String);
  void SetTargetLength(G4double);

  // -- everything split, with B field
  void makeSplitTrspTube();
  // -- decrepit: everything connected, but without B field
  void makeCombinedTrspTube();
  // -- accelerating field for atomic e-
  void makeAccel();
  // -- MCP and xtals
  void makeEndDetector();
  // -- target
  void makeTarget();
  // -- tracker
  void makeTracker();
  // -- endplate
  void makeEndplate();

private:
  int fVerbose;
  // -- Materials
  G4Material *fVac, *fAir, *fAl, *fC, *fBe, *fPb, *fSiO2,
    *fXeGas, *fAerog;

  // -- world
  G4Box*             fSolidWorld;
  G4LogicalVolume*   fLogicWorld;
  G4VPhysicalVolume* fPhysiWorld;

  // -- target
  G4Box*             fSolidTarget;
  G4LogicalVolume*   fLogicTarget;
  G4VPhysicalVolume* fPhysiTarget;
  G4Material*        fTargetMater;


  // -- endplates
  G4Box*             fSolidEndplate;
  G4LogicalVolume*   fLogicEndplate;
  G4VPhysicalVolume* fPhysiEndplate;
  G4Box*             fSolidEndplateF;
  G4LogicalVolume*   fLogicEndplateF;
  G4VPhysicalVolume* fPhysiEndplateF;

  // -- beampipe
  G4Tubs*            fSolidBeampipe;
  G4LogicalVolume*   fLogicBeampipe;
  G4VPhysicalVolume* fPhysiBeampipe;
  G4Material*        fBeampipeMater;


  // -- tracker
  G4Tubs*            fSolidTracker;
  G4LogicalVolume*   fLogicTracker;
  G4VPhysicalVolume* fPhysiTracker;

  // -- tracker chambers
  G4Box*             fSolidChamber;
  G4LogicalVolume**  fLogicChamber;
  G4VPhysicalVolume**fPhysiChamber;
  G4Material*        fChamberMater;

  G4Tubs*            fSolidShield;
  G4LogicalVolume*   fLogicShield;
  G4VPhysicalVolume* fPhysiShield;


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

  // -- End detector: MCP and xtals
  G4LogicalVolume*   fLogicMCP;
  G4VPhysicalVolume* fPhysiMCP;
  static const G4int fMacsXtalNum = 8;
  G4LogicalVolume**  fLogicXtals;
  G4VPhysicalVolume**fPhysiXtals;


  static G4ThreadLocal MagneticField* fpMagField;
  static G4ThreadLocal G4FieldManager* fpFieldMgr;
  static G4ThreadLocal G4FieldManager* fpFieldMgrE;

  DetectorMessenger* fDetectorMessenger;  // pointer to the Messenger
  G4UserLimits* fStepLimit;            // pointer to user step limits

  G4double fWorldLength;
  G4double fTargetLength, fBeampipeLength;
  G4double fBeampipeOuterRadius, fBeampipeInnerRadius, fTrkLength, fTrkOuterRadius, fTrkInnerRadius;
  G4int fNbOfChambers;
  G4double fChamberWidth;
  G4double fChamberSpacing;

  static const G4int fMacsTrkNum = 5;
  G4double fMacsTrkLength[fMacsTrkNum], fMacsTrkInnerRadius[fMacsTrkNum],
    fMacsTrkRadialThickness;

  G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

  //  G4Cache<ElectricFieldSetup*> fEmFieldSetup;
  G4Cache<F02ElectricFieldSetup*> fEmFieldSetup;

};


#endif
