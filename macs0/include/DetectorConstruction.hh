#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "MagneticField.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
class DetectorMessenger;
class G4GlobalMagFieldMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:

  DetectorConstruction();
  ~DetectorConstruction();

public:

  virtual G4VPhysicalVolume* Construct();
  void ConstructSDandField();
  G4VPhysicalVolume* macs1(); // MACS version 1, extracted from PRL,82,49

  const
  G4VPhysicalVolume* GetTracker() {return fPhysiTracker;};
  G4double GetTrackerFullLength() {return fTrkLength;};
  G4double GetTargetFullLength()  {return fTgtLength;};
  G4double GetWorldFullLength()   {return fWorldLength;};

  void SetTargetMaterial (G4String);
  void SetChamberMaterial(G4String);
  void SetMagField(G4double);

private:
  static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;


  G4Box*             fSolidWorld;    // pointer to the solid envelope
  G4LogicalVolume*   fLogicWorld;    // pointer to the logical envelope
  G4VPhysicalVolume* fPhysiWorld;    // pointer to the physical envelope

  G4Box*             fSolidTarget;   // pointer to the solid Target
  G4LogicalVolume*   fLogicTarget;   // pointer to the logical Target
  G4VPhysicalVolume* fPhysiTarget;   // pointer to the physical Target

  G4Tubs*            fSolidTracker;  // pointer to the solid Tracker
  G4LogicalVolume*   fLogicTracker;  // pointer to the logical Tracker
  G4VPhysicalVolume* fPhysiTracker;  // pointer to the physical Tracker

  G4Box*             fSolidChamber;  // pointer to the solid Chamber
  G4LogicalVolume**  fLogicChamber;  // pointer to the logical Chamber
  G4VPhysicalVolume**fPhysiChamber;  // pointer to the physical Chamber

  G4Material*         fTargetMater;  // pointer to the target  material
  G4Material*         fChamberMater; // pointer to the chamber material
  MagneticField* fPMagField;   // pointer to the magnetic field

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
};


#endif
