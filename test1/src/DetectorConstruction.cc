#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "ChamberParameterisation.hh"
#include "MagneticField.hh"
#include "TrackerSD.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

// ----------------------------------------------------------------------
DetectorConstruction::DetectorConstruction() :G4VUserDetectorConstruction(),
					      fSolidWorld(0), fLogicWorld(0), fPhysiWorld(0),
					      fSolidTarget(0), fLogicTarget(0), fPhysiTarget(0),
					      fSolidTracker(0), fLogicTracker(0), fPhysiTracker(0),
					      fSolidChamber(0), fPhysiChamber(0),
					      fTargetMater(0), fChamberMater(0), fPMagField(0), fDetectorMessenger(0),
					      fStepLimit(NULL),
					      fWorldLength(0.),
					      fTgtLength(0.), fTrkLength(0.), fTrkOuterRadius(0.),
					      fNbOfChambers(0), fChamberWidth(0.), fChamberSpacing(0.) {
  fPMagField = new MagneticField();
  fDetectorMessenger = new DetectorMessenger(this);
  fNbOfChambers = 2;
  fLogicChamber = new G4LogicalVolume*[fNbOfChambers];
}

// ----------------------------------------------------------------------
DetectorConstruction::~DetectorConstruction() {
  delete [] fLogicChamber;

  delete fPMagField;
  delete fDetectorMessenger;
}

// ----------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct() {

  // --------- Material definition ---------
  G4double a, z;
  G4double density, temperature, pressure;
  G4int nel;

  // -- Air
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);

  G4Material* Air = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);

  // -- Lead
  G4Material* Pb = new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

  // -- Aluminum
  G4Material* Al = new G4Material("Aluminum", z=13., a=26.98*g/mole, density=2.700*g/cm3);


  // -- Xenon gas
  G4Material* Xenon =  new G4Material("XenonGas", z=54., a=131.29*g/mole, density= 5.458*mg/cm3,
				      kStateGas, temperature= 293.15*kelvin, pressure= 1*atmosphere);

  // Print all the materials defined.
  //
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;


  fTgtLength  = 5.0*cm;                        // Full length of Target
  fTargetMater  = Pb;

  fNbOfChambers = 2;
  fChamberWidth = 2*cm;
  fChamberSpacing = 84*cm;
  fChamberMater = Xenon;

  fTrkLength = 200.*cm;
  fTrkOuterRadius = 100.*cm;

  fWorldLength= 1.5 * (fTgtLength + fTrkLength);

  G4double tgtHalfLength  = 0.5*fTgtLength;    // Half length of the Target
  G4double trkHalfLength = 0.5*fTrkLength;   // Half length of the Tracker

  // ------------------------------
  // -- World
  // ------------------------------
  G4double worldHalfLength = 0.5*fWorldLength;
  fSolidWorld = new G4Box("World", worldHalfLength, worldHalfLength, worldHalfLength);
  fLogicWorld = new G4LogicalVolume(fSolidWorld, Air, "World", 0, 0, 0);
  fPhysiWorld = new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0);

  // ------------------------------
  // -- Target
  // ------------------------------
  G4ThreeVector positionTarget = G4ThreeVector(0, 0, -(tgtHalfLength+trkHalfLength));
  fSolidTarget = new G4Box("target", tgtHalfLength, tgtHalfLength, tgtHalfLength);
  fLogicTarget = new G4LogicalVolume(fSolidTarget,fTargetMater,"Target",0,0,0);
  fPhysiTarget = new G4PVPlacement(0, positionTarget, fLogicTarget, "Target", fLogicWorld, false, 0);

  G4cout << "Target is " << fTgtLength/cm << "cm of " << fTargetMater->GetName() << G4endl;

  // ------------------------------
  // -- Tracker
  // ------------------------------
  G4ThreeVector positionTracker = G4ThreeVector(0,0,0);

  G4Tubs* trkTube = new G4Tubs("Tracker", 0, fTrkOuterRadius, trkHalfLength, 0.*deg, 360.*deg);
  fLogicTracker = new G4LogicalVolume(trkTube, Al, "Tracker", 0, 0, 0);
  fPhysiTracker = new G4PVPlacement(0, positionTracker, fLogicTracker, "Tracker", fLogicWorld, false, 0, true);

  // // Visualization attributes
  // G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  // G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));

  // fLogicWorld  ->SetVisAttributes(boxVisAtt);
  // fLogicTarget ->SetVisAttributes(boxVisAtt);
  // fLogicTracker->SetVisAttributes(boxVisAtt);


  // // Tracker segments

  // G4cout << "There are " << fNbOfChambers << " chambers in the tracker region. "
  //        << G4endl
  //        << "The chambers are " << fChamberWidth/cm << " cm of "
  //        << fChamberMater->GetName() << G4endl
  //        << "The distance between chamber is " << fChamberSpacing/cm << " cm"
  //        << G4endl;

  // G4double firstPosition = -trkHalfLength + fChamberSpacing;
  // G4double firstLength   = fTrkLength/10;
  // G4double lastLength    = fTrkLength;

  // G4double halfWidth = 0.5*fChamberWidth;
  // G4double rmaxFirst = 0.5 * firstLength;

  // G4double rmaxIncr = 0.0;
  // if( fNbOfChambers > 0 ){
  //   rmaxIncr =  0.5 * (lastLength-firstLength)/(fNbOfChambers-1);
  //   if (fChamberSpacing  < fChamberWidth) {
  //     G4Exception("B2aDetectorConstruction::DefineVolumes()",
  // 		  "InvalidSetup", FatalException,
  // 		  "Width>Spacing");
  //   }
  // }

  // for (G4int copyNo=0; copyNo<fNbOfChambers; copyNo++) {

  //   G4double Zposition = firstPosition + copyNo * fChamberSpacing;
  //   G4double rmax =  rmaxFirst + copyNo * rmaxIncr;

  //   G4Tubs* chamberS = new G4Tubs("Chamber_solid", 0, rmax, halfWidth, 0.*deg, 360.*deg);

  //   fLogicChamber[copyNo] =   new G4LogicalVolume(chamberS, fChamberMater, "Chamber_LV", 0, 0, 0);

  //   fLogicChamber[copyNo]->SetVisAttributes(chamberVisAtt);

  //   new G4PVPlacement(0,                            // no rotation
  // 		      G4ThreeVector(0,0,Zposition), // at (x,y,z)
  // 		      fLogicChamber[copyNo],        // its logical volume
  // 		      "Chamber_PV",                 // its name
  // 		      fLogicTracker,                    // its mother  volume
  // 		      false,                        // no boolean operations
  // 		      copyNo,                       // copy number
  // 		      true);              // checking overlaps

  // }

  // G4double maxStep = 0.5*fChamberWidth;
  // fStepLimit = new G4UserLimits(maxStep);
  // fLogicTracker->SetUserLimits(fStepLimit);





  // G4ThreeVector positionTracker = G4ThreeVector(0,0,0);

  // fSolidTracker = new G4Box("tracker", trkHalfLength, trkHalfLength, trkHalfLength);
  // fLogicTracker = new G4LogicalVolume(fSolidTracker , Air, "Tracker",0,0,0);
  // fPhysiTracker = new G4PVPlacement(0,              // no rotation
  // 				    positionTracker, // at (x,y,z)
  // 				    fLogicTracker,    // its logical volume
  // 				    "Tracker",       // its name
  // 				    fLogicWorld,      // its mother  volume
  // 				    false,           // no boolean operations
  // 				    0);              // copy number

  // // ----------------------------------------------------------------------
  // // Tracker segments with parametrized volumes
  // // ----------------------------------------------------------------------
  // fSolidChamber = new G4Box("chamber", 100*cm, 100*cm, 10*cm);
  // fLogicChamber = new G4LogicalVolume(fSolidChamber, fChamberMater,"Chamber",0,0,0);

  // G4double firstPosition = -trkHalfLength + 0.5*fChamberWidth;
  // G4double firstLength = fTrkLength/10;
  // G4double lastLength  = fTrkLength;

  // G4VPVParameterisation* chamberParam = new ChamberParameterisation(
  // 								    fNbOfChambers,          // NoChambers
  // 								    firstPosition,         // Z of center of first
  // 								    fChamberSpacing,        // Z spacing of centers
  // 								    fChamberWidth,          // Width Chamber
  // 								    firstLength,           // lengthInitial
  // 								    lastLength);           // lengthFinal

  // fPhysiChamber = new G4PVParameterised(
  // 					"Chamber",       // their name
  // 					fLogicChamber,    // their logical volume
  // 					fLogicTracker,    // Mother logical volume
  // 					kZAxis,          // Are placed along this axis
  // 					fNbOfChambers,    // Number of chambers
  // 					chamberParam);   // The parametrisation

  // G4cout << "There are " << fNbOfChambers << " chambers in the tracker region. "
  //        << "The chambers are " << fChamberWidth/mm << " mm of "
  //        << fChamberMater->GetName() << "\n The distance between chamber is "
  //        << fChamberSpacing/cm << " cm" << G4endl;

  // //------------------------------------------------
  // // Sensitive detectors
  // //------------------------------------------------

  // G4SDManager* SDman = G4SDManager::GetSDMpointer();

  // G4String trackerChamberSDname = "/TrackerChamberSD";
  // TrackerSD* aTrackerSD = new TrackerSD( trackerChamberSDname );
  // SDman->AddNewDetector( aTrackerSD );
  // fLogicChamber->SetSensitiveDetector( aTrackerSD );

  // //--------- Visualization attributes -------------------------------

  // G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  // fLogicWorld  ->SetVisAttributes(BoxVisAtt);
  // fLogicTarget ->SetVisAttributes(BoxVisAtt);
  // fLogicTracker->SetVisAttributes(BoxVisAtt);

  // G4VisAttributes* ChamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  // fLogicChamber->SetVisAttributes(ChamberVisAtt);

  // G4double maxStep = 0.5*fChamberWidth;
  // fLogicTracker->SetUserLimits(new G4UserLimits(maxStep));

  return fPhysiWorld;
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
