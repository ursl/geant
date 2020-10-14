#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "ChamberParameterisation.hh"
#include "MagneticField.hh"
#include "TrackerSD.hh"

#include "TTree.h"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


G4ThreadLocal G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = 0;

// ----------------------------------------------------------------------
DetectorConstruction::DetectorConstruction() :G4VUserDetectorConstruction(),
					      fSolidWorld(0), fLogicWorld(0), fPhysiWorld(0),
					      fSolidTarget(0), fLogicTarget(0), fPhysiTarget(0),
					      fSolidTracker(0), fLogicTracker(0), fPhysiTracker(0),
					      fSolidChamber(0),
					      fTargetMater(0), fChamberMater(0), fPMagField(0), fDetectorMessenger(0),
					      fStepLimit(NULL),
					      fWorldLength(0.),
					      fTgtLength(0.), fTrkLength(0.), fTrkOuterRadius(0.), fTrkInnerRadius(0.),
					      fNbOfChambers(0), fChamberWidth(0.), fChamberSpacing(0.), fCheckOverlaps(true)
{
  fPMagField = new MagneticField();
  fDetectorMessenger = new DetectorMessenger(this);
  fNbOfChambers = 2;
  fLogicChamber = new G4LogicalVolume*[fMacsTrkNum];
  fPhysiChamber = new G4VPhysicalVolume*[fMacsTrkNum];
}

// ----------------------------------------------------------------------
DetectorConstruction::~DetectorConstruction() {
  delete []fLogicChamber;
  delete []fPhysiChamber;

  delete fPMagField;
  delete fDetectorMessenger;
}

// ----------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct() {
  return macs1();
}

// ----------------------------------------------------------------------
void DetectorConstruction::SetTargetMaterial(G4String materialName) {
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);
  if (pttoMaterial)
    {fTargetMater = pttoMaterial;
      fLogicTarget->SetMaterial(pttoMaterial);
      G4cout << "\n----> The target is " << fTgtLength/cm << " cm of "
             << materialName << G4endl;
    }
}

// ----------------------------------------------------------------------
void DetectorConstruction::SetChamberMaterial(G4String materialName) {
  G4NistManager* nistManager = G4NistManager::Instance();

  G4Material* pttoMaterial =
              nistManager->FindOrBuildMaterial(materialName);

  if (fChamberMater != pttoMaterial) {
     if ( pttoMaterial ) {
        fChamberMater = pttoMaterial;
        for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {
            if (fLogicChamber[copyNo]) fLogicChamber[copyNo]->SetMaterial(fChamberMater);
        }
        G4cout
          << G4endl
          << "----> The chambers are made of " << materialName << G4endl;
     } else {
        G4cout
          << G4endl
          << "-->  WARNING from SetChamberMaterial : "
          << materialName << " not found" << G4endl;
     }
  }
}

// ----------------------------------------------------------------------
void DetectorConstruction::SetMagField(G4double fieldValue) {
  fPMagField->SetFieldValue(fieldValue);
}



// ----------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::macs1() {

  G4double a, z;
  G4double density, temperature, pressure;
  G4int nel, natoms;

  // -- Air
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);

  G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);

  // -- SiO2
  G4Element* Si    = new G4Element("Silicon", "Si", z=14., a = 28.09*g/mole);
  G4Material* SiO2 = new G4Material("SiO2", density = 2.5*g/cm3, nel = 2);
  SiO2->AddElement (Si, natoms=1);
  SiO2->AddElement (O, natoms=2);

  // -- Aluminum
  G4Material* Al = new G4Material("Aluminum", z=13., a=26.98*g/mole, density=2.700*g/cm3);

  // -- Lead
  G4Material* Pb = new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

  // -- Xenon gas
  G4Material* Xenon =  new G4Material("XenonGas", z=54., a=131.29*g/mole, density= 5.458*mg/cm3,
				      kStateGas, temperature= 293.15*kelvin, pressure= 1*atmosphere);

  // -- placeholder
  fChamberMater = Xenon;

  // -- MACS detector dimensions
  fMacsTrkRadialThickness = 1.0*cm; // placeholder
  for (int i = 0; i < fMacsTrkNum; ++i) {
    fMacsTrkLength[i] = 38.0*cm + i*10.0*cm;
    fMacsTrkInnerRadius[i] = 6.0*cm + i*5.0*cm;
  }

  // -- target thickness: d = 8mg/cm2 -> l = d/rho = (8mg/cm) / (2.5e3mg/cm3) = 3.2e-3cm
  fTgtLength  = 3.2e-3*cm;                        // Full length of Target
  fTgtLength  = 3.2*cm;                        // placeholder
  fTargetMater  = SiO2;
  fTargetMater  = Pb;                          // placeholder
  fTgtLength  = 0.1*cm;                        // placeholder

  fWorldLength= 200.0*cm;

  fTrkOuterRadius = 100.0*cm;
  G4double tgtHalfLength  = 0.5*fTgtLength;    // Half length of the Target
  G4double trkHalfLength  = 0.5*100.*cm;

  // ------------------------------
  // -- World
  // ------------------------------
  G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  G4double worldHalfLength = 0.5*fWorldLength;
  fSolidWorld = new G4Box("World", worldHalfLength, worldHalfLength, worldHalfLength);
  fLogicWorld = new G4LogicalVolume(fSolidWorld, Air, "World", 0, 0, 0);
  fPhysiWorld = new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0);

  // ------------------------------
  // -- Target
  // ------------------------------
  G4ThreeVector positionTarget = G4ThreeVector(0., 0., 0.);
  fSolidTarget = new G4Box("Target", 5.0*cm, 5.0*cm, tgtHalfLength);
  fLogicTarget = new G4LogicalVolume(fSolidTarget,fTargetMater,"Target", 0, 0, 0);
  fPhysiTarget = new G4PVPlacement(0, positionTarget, fLogicTarget, "Target", fLogicWorld, false, 0);
  fLogicTarget->SetVisAttributes(boxVisAtt);

  G4cout << "Target is " << fTgtLength/cm << "cm of " << fTargetMater->GetName() << G4endl;

  // ------------------------------
  // -- Tracker
  // ------------------------------
  G4ThreeVector positionTracker = G4ThreeVector(0,0,0);

  fSolidTracker = new G4Tubs("Tracker", 0, fTrkOuterRadius, trkHalfLength, 0.*deg, 360.*deg);
  fLogicTracker = new G4LogicalVolume(fSolidTracker, Al, "Tracker", 0, 0, 0);
  fPhysiTracker = new G4PVPlacement(0, positionTracker, fLogicTracker, "Tracker", fLogicWorld, false, 0, true);
  fLogicTracker->SetVisAttributes(boxVisAtt);

  G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  for (G4int itrk=0; itrk < fMacsTrkNum; itrk++) {

    G4Tubs* chamberS = new G4Tubs("Chamber_solid", fMacsTrkInnerRadius[itrk],
				  fMacsTrkInnerRadius[itrk] + fMacsTrkRadialThickness,
				  0.5*fMacsTrkLength[itrk],
				  0.*deg, 360.*deg);

    fLogicChamber[itrk] = new G4LogicalVolume(chamberS, fChamberMater, Form("Chamber_LV"), 0, 0, 0);
    fLogicChamber[itrk]->SetVisAttributes(chamberVisAtt);

    fPhysiChamber[itrk] = new G4PVPlacement(0,                            // no rotation
					    G4ThreeVector(0., 0., 0.),    // at (x,y,z)
					    fLogicChamber[itrk],          // its logical volume
					    Form("Chamber_PV"),           // its name
					    fLogicTracker,                // its mother  volume
					    false,                        // no boolean operations
					    itrk,                         // copy number
					    fCheckOverlaps);              // checking overlaps

  }




  return fPhysiWorld;
}

// ----------------------------------------------------------------------
void DetectorConstruction::ConstructSDandField() {

  // Sensitive detectors
  TrackerSD* aTrackerSD = new TrackerSD("muamu/TrackerChamberSD", "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name of "Chamber_LV".
  SetSensitiveDetector("Chamber_LV", aTrackerSD, true);

  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue = G4ThreeVector();
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);

}
