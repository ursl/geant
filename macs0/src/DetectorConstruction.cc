#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "ChamberParameterisation.hh"
#include "MagneticField.hh"
#include "ElectricFieldSetup.hh"
#include "TrackerSD.hh"
#include "MCPSD.hh"

#include "TTree.h"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"
#include "G4Torus.hh"
#include "G4UnionSolid.hh"
#include "G4MultiUnion.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4SDManager.hh"
#include "G4NistManager.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"


G4ThreadLocal MagneticField* DetectorConstruction::fpMagField = 0;
G4ThreadLocal G4FieldManager* DetectorConstruction::fpFieldMgr = 0;


// ----------------------------------------------------------------------
DetectorConstruction::DetectorConstruction() :G4VUserDetectorConstruction(),
					      fSolidWorld(0), fLogicWorld(0), fPhysiWorld(0),
					      fSolidTarget(0), fLogicTarget(0), fPhysiTarget(0), fTargetMater(0),
					      fSolidTracker(0), fLogicTracker(0), fPhysiTracker(0),
					      fSolidChamber(0), fChamberMater(0),
					      fSolidTrsp(0), fLogicTrsp(0), fPhysiTrsp(0), fTrspMater(0),
					      fTrspOuterRadius(10.0), fTrspLength1(50.), fTrspLength2(50.),
					      fMagneticLogical(0),
					      fLogicAccel(0), fPhysiAccel(0),
					      fLogicMCP(0), fPhysiMCP(0),
					      fDetectorMessenger(0),
					      fStepLimit(NULL),
					      fWorldLength(0.),
					      fTgtLength(0.), fTrkLength(0.), fTrkOuterRadius(0.), fTrkInnerRadius(0.),
					      fNbOfChambers(0), fChamberWidth(0.), fChamberSpacing(0.), fCheckOverlaps(true)
{
  fDetectorMessenger = new DetectorMessenger(this);
  fNbOfChambers = 2;
  fLogicChamber = new G4LogicalVolume*[fMacsTrkNum];
  fPhysiChamber = new G4VPhysicalVolume*[fMacsTrkNum];
  fLogicXtals   = new G4LogicalVolume*[fMacsXtalNum];
  fPhysiXtals   = new G4VPhysicalVolume*[fMacsXtalNum];
}

// ----------------------------------------------------------------------
DetectorConstruction::~DetectorConstruction() {
  delete []fLogicChamber;
  delete []fPhysiChamber;

  delete fpMagField;
  delete fDetectorMessenger;
}

// ----------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct() {
  defineMaterials();
  return macs0();
}



// ----------------------------------------------------------------------
void DetectorConstruction::defineMaterials() {
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 3.e-18*pascal;
  // -- vacuum defined as intergalactic
  fVac =  new G4Material("vac", atomicNumber,  massOfMole, density, kStateGas, temperature, pressure);


  G4double a, z, fractionmass;
  G4int nel, natoms, ncomp;

  // -- Air
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole);

  fAir = new G4Material("Air", density= 1.29*mg/cm3, nel=2);
  fAir->AddElement(N, 70*perCent);
  fAir->AddElement(O, 30*perCent);

  // -- SiO2
  G4Element* Si    = new G4Element("Silicon", "Si", z=14., a = 28.09*g/mole);
  fSiO2 = new G4Material("SiO2", density = 2.5*g/cm3, nel = 2);
  fSiO2->AddElement (Si, natoms=1);
  fSiO2->AddElement (O, natoms=2);


  // -- Aerogel (from https://geant4.web.cern.ch/.../geometry/training/D1-Materials.pdf)
  G4Element* elH = new G4Element("Hydrogen", "H", z=1., 1.01*g/mole);
  G4Element* elC = new G4Element("Carbon", "C", 6., 12.011*g/mole);
  G4Element* elO = new G4Element("Oxygen", "O", z=8., 16.00*g/mole);
  G4Material* H2O = new G4Material("Water", density=1.000*g/cm3, ncomp=2);
  H2O->AddElement(elH, natoms=2);
  H2O->AddElement(elO, natoms=1);

  density = 0.200*g/cm3;
  fAerog = new G4Material("Aerogel", density, ncomp=3);
  fAerog->AddMaterial(fSiO2,fractionmass=62.5*perCent);
  fAerog->AddMaterial(H2O ,fractionmass=37.4*perCent);
  fAerog->AddElement (elC ,fractionmass= 0.1*perCent);

  // -- Aluminum
  fAl = new G4Material("Aluminum", z=13., a=26.98*g/mole, density=2.700*g/cm3);

  // -- Lead
  fPb = new G4Material("Lead", z=82., a= 207.19*g/mole, density= 11.35*g/cm3);

  // -- Xenon gas
  fXeGas =  new G4Material("XenonGas", z=54., a=131.29*g/mole, density= 5.458*mg/cm3,
			   kStateGas, temperature= 293.15*kelvin, pressure= 1*atmosphere);

}

// ----------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::macs0() {

  // -- placeholder
  fChamberMater = fXeGas;

  fTargetMater  = fAerog;                          // placeholder
  fTgtLength  = 5*cm;                        // placeholder
  fTrkLength  = 100.*cm;

  // -- MACS detector dimensions
  fChamberWidth = fMacsTrkRadialThickness = 1.0*cm; // placeholder
  for (int i = 0; i < fMacsTrkNum; ++i) {
    fMacsTrkLength[i] = fTrkLength; // 38.0*cm + (fMacsTrkNum-1)*10.0*cm;
    fMacsTrkInnerRadius[i] = 6.0*cm + i*5.0*cm;
  }


  fWorldLength= 420.0*cm;
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
  G4cout << "==========> Computed tolerance = "
         << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
         << " mm" << G4endl;

  fTrkOuterRadius = 100.0*cm;
  G4double tgtHalfLength  = 0.5*fTgtLength;    // Half length of the Target
  G4double trkHalfLength  = 0.5*fTrkLength;

  // ------------------------------
  // -- World
  // ------------------------------
  G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  G4double worldHalfLength = 0.5*fWorldLength;
  fSolidWorld = new G4Box("World", worldHalfLength, worldHalfLength, worldHalfLength);
  fLogicWorld = new G4LogicalVolume(fSolidWorld, fVac, "World", 0, 0, 0);
  fPhysiWorld = new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0, true);

  // ------------------------------
  // -- Target
  // ------------------------------
  //  G4ThreeVector positionTarget = G4ThreeVector(0., 0., -trkHalfLength-tgtHalfLength);
  G4ThreeVector positionTarget = G4ThreeVector(0., 0., -trkHalfLength);
  fSolidTarget = new G4Box("Target", 2*cm, 2*cm, tgtHalfLength);
  fLogicTarget = new G4LogicalVolume(fSolidTarget,fTargetMater,"Target", 0, 0, 0);
  fPhysiTarget = new G4PVPlacement(0, positionTarget, fLogicTarget, "Target", fLogicWorld, false, 0, true);
  fLogicTarget->SetVisAttributes(boxVisAtt);

  G4cout << "Target is " << fTgtLength/cm << "cm of " << fTargetMater->GetName() << G4endl;

  // ------------------------------
  // -- Tracker
  // ------------------------------
  G4ThreeVector positionTracker = G4ThreeVector(0,0,0);

  fSolidTracker = new G4Tubs("Tracker", 0, fTrkOuterRadius, trkHalfLength, 0.*deg, 360.*deg);
  fLogicTracker = new G4LogicalVolume(fSolidTracker, fVac, "Tracker", 0, 0, 0);
  fPhysiTracker = new G4PVPlacement(0, positionTracker, fLogicTracker, "Tracker", fLogicWorld, false, 0, true);
  fLogicTracker->SetVisAttributes(boxVisAtt);

  G4VisAttributes* chamberVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  for (G4int itrk=0; itrk < fMacsTrkNum; itrk++) {

    G4Tubs* chamberS = new G4Tubs("Chamber_solid", fMacsTrkInnerRadius[itrk],
				  fMacsTrkInnerRadius[itrk] + fMacsTrkRadialThickness,
				  0.5*fMacsTrkLength[itrk],
				  0.*deg, 360.*deg);

    fLogicChamber[itrk] = new G4LogicalVolume(chamberS, fChamberMater, Form("Chamber_LV"), 0, 0, 0);
    double iscale = 0.3 + itrk*0.1;
    fLogicChamber[itrk]->SetVisAttributes(G4Colour(0.8, iscale, iscale));

    fPhysiChamber[itrk] = new G4PVPlacement(0,                            // no rotation
					    G4ThreeVector(0., 0., 0.),    // at (x,y,z)
					    fLogicChamber[itrk],          // its logical volume
					    Form("Chamber_PV"),           // its name
					    fLogicTracker,                // its mother  volume
					    false,                        // no boolean operations
					    itrk,                         // copy number
					    fCheckOverlaps);              // checking overlaps

  }


  G4double maxStep = 0.5*fChamberWidth;
  fStepLimit = new G4UserLimits(maxStep);
  fLogicTracker->SetUserLimits(fStepLimit);


  // ------------------------------
  // -- Transport tube
  // ------------------------------
  makeSplitTrspTube();
  makeAccel();
  makeEndDetector();

  return fPhysiWorld;
}

// ----------------------------------------------------------------------
void DetectorConstruction::SetMagField(G4double fieldValue) {
  G4cout << "===============> DetectorConstruction::SetMagField> SetMagField(" << fieldValue << G4endl;
  fpMagField->SetField(fieldValue);
}


// ----------------------------------------------------------------------
void DetectorConstruction::ConstructSDandField() {

  // Sensitive detectors
  TrackerSD* aTrackerSD = new TrackerSD("macs0/TrackerChamberSD", "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  SetSensitiveDetector("Chamber_LV", aTrackerSD, true);

  MCPSD* aMCPSD = new MCPSD("macs0/MCPSD", "MCPHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aMCPSD);
  SetSensitiveDetector("lMCP", aMCPSD, true);


  fpMagField = new MagneticField();
  fpFieldMgr = new G4FieldManager();
  fpFieldMgr->SetDetectorField(fpMagField);
  fpFieldMgr->CreateChordFinder(fpMagField);
  fMagneticLogical->SetFieldManager(fpFieldMgr, true);

  ElectricFieldSetup* fieldSetup = new ElectricFieldSetup();
  fEmFieldSetup.Put(fieldSetup);
  fLogicAccel->SetFieldManager(fpFieldMgr, true);


  // -- these two lines lead to BREAKs
  //  G4AutoDelete::Register(fpMagField);
  //  G4AutoDelete::Register(fpFieldMgr);
}

// ----------------------------------------------------------------------
void DetectorConstruction::makeAccel() {

  G4double orAcc = 26.0*cm;
  G4double lAcc =   2.0*cm;
  auto p0 = new G4Tubs("Accel", orAcc-0.1*cm, orAcc, lAcc, 0.*deg, 360.*deg);
  fLogicAccel = new G4LogicalVolume(p0, fVac, "lAccel");

  G4RotationMatrix* fieldRot = new G4RotationMatrix();
  //  fieldRot->rotateX(90.*deg);
  fPhysiAccel =  new G4PVPlacement(fieldRot,
				   G4ThreeVector(0., 0., 0.5*(fTrkLength+2.*lAcc)+0.5*cm),
				   fLogicAccel, "pAccel", fLogicWorld, false, 0, fCheckOverlaps);

  G4UserLimits* userLimits = new G4UserLimits(0.1*cm);
  fLogicAccel->SetUserLimits(userLimits);

  G4VisAttributes *pVA  = new G4VisAttributes;
  pVA->SetColour(G4Colour(0.0, 0.8, 0.8));
  pVA->SetForceSolid(true);
  fLogicAccel->SetVisAttributes(pVA);

}

// ----------------------------------------------------------------------
void DetectorConstruction::makeEndDetector() {

  G4double orMCP = 6.0*cm;
  G4double lMCP =   0.1*cm;
  auto p0 = new G4Tubs("MCP", 0., orMCP, lMCP, 0.*deg, 360.*deg);
  fLogicMCP = new G4LogicalVolume(p0, fAl, "lMCP");

  G4RotationMatrix* fieldRot = new G4RotationMatrix();
  fieldRot->rotateY(90.*deg);
  fPhysiMCP =  new G4PVPlacement(fieldRot,
				 G4ThreeVector(-(fTrspLength2 + 40.0*cm + 2.0*cm), 0., fTrkLength + fTrspLength1 - 5.0*cm),
				 fLogicMCP, "pMCP", fLogicWorld, false, 0, fCheckOverlaps);

  // set step limit in tube with magnetic field
  G4UserLimits* userLimits = new G4UserLimits(0.1*cm);
  fLogicMCP->SetUserLimits(userLimits);

  G4VisAttributes *pVA  = new G4VisAttributes;
  pVA->SetColour(G4Colour(0.6, 0.4, 0.2));
  pVA->SetForceWireframe(true);
  fLogicMCP->SetVisAttributes(pVA);


}


// ----------------------------------------------------------------------
void DetectorConstruction::makeSplitTrspTube() {

  fTrspOuterRadius = 6.*cm;
  fTrspLength1 = fTrspLength2 = 60.*cm;
  G4double pRtor = 40.*cm;
  fTrspMater  = fAl;
  G4double trkHalfLength  = 0.5*fTrkLength;

  G4ThreeVector positionTrsp = G4ThreeVector(0, 0, 2.*trkHalfLength - 0.25*fTrspLength1);
  G4Tubs  *part1 = new G4Tubs("Trsp1", fTrspOuterRadius-0.1*cm, fTrspOuterRadius, 0.5*fTrspLength1, 0.*deg, 360.*deg);
  G4Tubs  *part3 = new G4Tubs("Trsp3", fTrspOuterRadius-0.1*cm, fTrspOuterRadius, 0.5*fTrspLength2, 0.*deg, 360.*deg);


  // -- trafo for part3
  G4RotationMatrix* rot3 = new G4RotationMatrix;
  rot3->rotateY(M_PI/2.*rad);
  G4ThreeVector zTrans3(pRtor + 0.5*fTrspLength2, 0, -0.5*fTrspLength1 - pRtor);
  G4RotationMatrix invRot3 = rot3->invert();
  G4Transform3D transform3(invRot3, -zTrans3);

  G4Transform3D unity;

  fSolidTrsp = new G4MultiUnion("allTrsp");
  fSolidTrsp->AddNode(*part1, unity);
  fSolidTrsp->AddNode(*part3, transform3);
  fSolidTrsp->Voxelize();

  fLogicTrsp = new G4LogicalVolume(fSolidTrsp, fAl, "Trsp", 0, 0, 0);
  fPhysiTrsp = new G4PVPlacement(0, positionTrsp, fLogicTrsp, "Trsp", fLogicWorld, false, 0, true);

  G4VisAttributes *pVA  = new G4VisAttributes;
  pVA->SetColour(G4Colour(0.8, 0.8, 0.8));
  pVA->SetForceSolid(true);
  fLogicTrsp->SetVisAttributes(pVA);

  // -- now add tube with magnetic field, just as in B5:
  G4double fdiam = 0.36*m;
  G4double fheight = 0.31*m;
  auto magneticSolid = new G4Tubs("magneticTubs", 0., fdiam, fheight, 0., 360.*deg);
  fMagneticLogical = new G4LogicalVolume(magneticSolid, fVac, "magneticLogical");

  G4RotationMatrix* fieldRot = new G4RotationMatrix();
  fieldRot->rotateX(90.*deg);

  G4ThreeVector zTransField(0, 0, -fTrspLength1 - fTrkLength);
  G4RotationMatrix invFieldRot = fieldRot->invert();
  G4Transform3D transformField(invFieldRot, -zTransField);


  new G4PVPlacement(transformField, fMagneticLogical, "magneticPhysical", fLogicWorld, false, 0, true);

  // set step limit in tube with magnetic field
  G4UserLimits* userLimits = new G4UserLimits(0.1*cm);
  fMagneticLogical->SetUserLimits(userLimits);

}


// ----------------------------------------------------------------------
void DetectorConstruction::makeCombinedTrspTube() {

  fTrspOuterRadius = 6.*cm;
  fTrspLength1 = fTrspLength2 = 60.*cm;
  G4double pRtor = 40.*cm;
  fTrspMater  = fAl;
  G4double trkHalfLength  = 0.5*fTrkLength;

  G4ThreeVector positionTrsp = G4ThreeVector(0, 0, 2.*trkHalfLength - 0.25*fTrspLength1);
  G4Tubs  *part1 = new G4Tubs("Trsp1", fTrspOuterRadius-0.1*cm, fTrspOuterRadius, 0.5*fTrspLength1, 0.*deg, 360.*deg);
  G4Torus *part2 = new G4Torus("Trsp2", fTrspOuterRadius-0.1*cm, fTrspOuterRadius, pRtor, 0.*degree, 90.*degree);
  fMagneticLogical = new G4LogicalVolume(part2, fVac, "magneticLogical");
  G4Tubs  *part3 = new G4Tubs("Trsp3", fTrspOuterRadius-0.1*cm, fTrspOuterRadius, 0.5*fTrspLength2, 0.*deg, 360.*deg);

  // -- trafo for part2
  G4RotationMatrix* rot2 = new G4RotationMatrix;
  rot2->rotateX(-M_PI/2.*rad);
  G4ThreeVector zTrans2(pRtor, 0, -0.5*fTrspLength1);
  G4RotationMatrix invRot2 = rot2->invert();
  G4Transform3D transform2(invRot2, -zTrans2);

  // -- trafo for part3
  G4RotationMatrix* rot3 = new G4RotationMatrix;
  rot3->rotateY(M_PI/2.*rad);
  G4ThreeVector zTrans3(pRtor + 0.5*fTrspLength2, 0, -0.5*fTrspLength1 - pRtor);
  G4RotationMatrix invRot3 = rot3->invert();
  G4Transform3D transform3(invRot3, -zTrans3);

  G4Transform3D unity;

  fSolidTrsp = new G4MultiUnion("allTrsp");
  fSolidTrsp->AddNode(*part1, unity);
  fSolidTrsp->AddNode(*part2, transform2);
  fSolidTrsp->AddNode(*part3, transform3);
  fSolidTrsp->Voxelize();

  fLogicTrsp = new G4LogicalVolume(fSolidTrsp, fAl, "Trsp", 0, 0, 0);
  fPhysiTrsp = new G4PVPlacement(0, positionTrsp, fLogicTrsp, "Trsp", fLogicWorld, false, 0, true);

  G4VisAttributes *pVA  = new G4VisAttributes;
  pVA->SetColour(G4Colour(0.8, 0.8, 0.8));
  pVA->SetForceSolid(true);
  fLogicTrsp->SetVisAttributes(pVA);
}
