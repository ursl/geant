#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <math.h>

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "ChamberParameterisation.hh"
#include "MagneticField.hh"
#include "ElectricFieldSetup.hh"
#include "RegionInformation.hh"
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

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"

#include "G4GeometryTolerance.hh"
#include "G4UserLimits.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

using namespace std;


G4ThreadLocal MagneticField* DetectorConstruction::fpMagField = 0;
G4ThreadLocal G4FieldManager* DetectorConstruction::fpFieldMgr = 0;


// ----------------------------------------------------------------------
DetectorConstruction::DetectorConstruction() :G4VUserDetectorConstruction(),
					      fSolidWorld(0), fLogicWorld(0), fPhysiWorld(0),
					      fSolidTarget(0), fLogicTarget(0), fPhysiTarget(0), fTargetMater(0),
					      fSolidBeampipe(0), fLogicBeampipe(0), fPhysiBeampipe(0),
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
					      fTargetLength(0.), fBeampipeLength(0.),
					      fBeampipeInnerRadius(5.),fBeampipeOuterRadius(5.02),
					      fTrkLength(0.), fTrkInnerRadius(5.), fTrkOuterRadius(100.),
					      fNbOfChambers(0), fChamberWidth(0.), fChamberSpacing(0.),
					      fCheckOverlaps(true)
{
  G4cout << "==DetectorConstruction> c'tor called"
	 << G4endl;
  fDetectorMessenger = new DetectorMessenger(this);
  fNbOfChambers = 2;
  fLogicChamber = new G4LogicalVolume*[fMacsTrkNum];
  fPhysiChamber = new G4VPhysicalVolume*[fMacsTrkNum];
  fLogicXtals   = new G4LogicalVolume*[fMacsXtalNum];
  fPhysiXtals   = new G4VPhysicalVolume*[fMacsXtalNum];

  defineMaterials();

}

// ----------------------------------------------------------------------
DetectorConstruction::~DetectorConstruction() {
  delete []fLogicChamber;
  delete []fPhysiChamber;

  delete fpMagField;
  delete fDetectorMessenger;
}

// ----------------------------------------------------------------------
void DetectorConstruction::parseCmd(G4String filename) {

  ifstream INS;
  string sline;
  INS.open(filename);
  string tmat("/det/setTargetMaterial ");
  string tlen("/det/setTargetLength ");
  while (getline(INS, sline)) {
    if (0 == sline.find("#")) continue;
    // -- search for target material
    if (string::npos != sline.find(tmat)) {
      size_t start_pos = sline.find(tmat) + tmat.size();
      cout << "Found tmat at position " << start_pos << endl;
      string material = sline.substr(start_pos);
      cout << "material ->" << material << "<-" << endl;
      SetTargetMaterial(material);
    }
    // -- search for target length
    if (string::npos != sline.find(tlen)) {
      size_t start_pos = sline.find(tlen) + tlen.size();
      cout << "Found tlen at position " << start_pos << endl;
      string length = sline.substr(start_pos);
      cout << "length ->" << length << "<-" << endl;

      start_pos = length.find(" ");
      string unit = length.substr(start_pos + 1);
      //      SetTargetLength(lenmm);
      cout << "unit ->" << unit << "<-" << endl;

      string number = length.substr(0, length.find(unit)-1);
      cout << "number ->" <<  number << "<-" << endl;
      double s(1.);
      if ("mm" == unit) {
	s = 1.;
      } else if ("cm" == unit) {
	s = 10.;
      } else if ("nm" == unit) {
	s = 1.e-6;
      } else if ("um" == unit) {
	s = 1.e-3;
      }
      G4double unitnumber = stod(number)*s;
      SetTargetLength(unitnumber);
      cout << "setTargetLength(" << unitnumber << ")" << endl;
    }
  }

  INS.close();
}

// ----------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct() {
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  G4cout << "==DetectorConstruction::Construct>  building detector"
	 << G4endl;
  return macs0();
}


// ----------------------------------------------------------------------
void DetectorConstruction::SetTargetLength(G4double value) {

  fTargetLength = value;
  G4cout << "==DetectorConstruction::SetTargetMaterial> (re)setting target length to " << fTargetLength
	 << G4endl;
}


// ----------------------------------------------------------------------
void DetectorConstruction::SetTargetMaterial(G4String materialName) {

  // -- if custom is set, material defined in this class is used!
  bool custom(false);

  if (materialName == "Aerogel") {
    fTargetMater = fAerog;
    custom = true;
  } else if (materialName == "SiO2") {
    fTargetMater = fSiO2;
    custom = true;
  } else if (materialName == "Cfoil") {
    fTargetMater = fC;
    custom = true;
  }

  if (custom) return;

  G4NistManager *nistManager = G4NistManager::Instance();
  G4Material *pMaterial = nistManager->FindOrBuildMaterial(materialName);

  if (pMaterial) {
    fTargetMater = pMaterial;
    G4cout << "==DetectorConstruction::SetTargetMaterial> The target is made of " << materialName << G4endl;
  } else {
    G4cout << "==DetectorConstruction::SetTargetMaterial> XXX  WARNING from SetTargetMaterial : "
	   << materialName << " not found XXX" << G4endl;
  }
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
  fBe = new G4Material("Beryllium", z=4, a=9.012182*g/mole, density=1.848*g/cm3);

  // -- Carbon (fiber? foil??)
  density = 0.145*g/cm3;
  fC = new G4Material("CarbonFoil", density, nel=1);
  fC->AddElement(elC, 1);

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
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // -- placeholder
  fChamberMater = fXeGas;

  fBeampipeMater  = fBe;
  // -- provide default (can be overridden with mac)
  if (0) {
    fTargetMater  = fAerog;
    fTargetLength  = 1*cm;

    fTargetMater  = fAl;
    fTargetLength  = 0.1*mm; // 10um

    // -- EPJ, C80, 804:
    fTargetMater  = fC;
    fTargetLength  = 15*nm; // 15nm
  }

  //
  fTrkLength  = 100.*cm;
  fBeampipeLength = 120.*cm;

  // -- MACS detector dimensions
  fChamberWidth = fMacsTrkRadialThickness = 1.0*cm; // placeholder
  for (int i = 0; i < fMacsTrkNum; ++i) {
    fMacsTrkLength[i] = fTrkLength; // 38.0*cm + (fMacsTrkNum-1)*10.0*cm;
    fMacsTrkInnerRadius[i] = 6.0*cm + i*5.0*cm;
  }


  fWorldLength= 420.0*cm;
  // G4GeometryManager::GetInstance()->SetWorldMaximumExtent(fWorldLength);
  // G4cout << "==========> Computed tolerance = "
  //        << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
  //        << " mm" << G4endl;

  fBeampipeInnerRadius = 3.0*cm;
  fBeampipeOuterRadius = 3.1*cm;


  fTrkInnerRadius = 5.0*cm;
  fTrkOuterRadius = 100.0*cm;
  G4double trkHalfLength  = 0.5*fTrkLength;
  G4double beampipeHalfLength  = 0.5*fBeampipeLength;

  // ------------------------------
  // -- World
  // ------------------------------
  G4VisAttributes* boxVisAtt= new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  G4double worldHalfLength = 0.5*fWorldLength;
  fSolidWorld = new G4Box("World", worldHalfLength, worldHalfLength, worldHalfLength);
  fLogicWorld = new G4LogicalVolume(fSolidWorld, fVac, "World", 0, 0, 0);
  fPhysiWorld = new G4PVPlacement(0, G4ThreeVector(), fLogicWorld, "World", 0, false, 0, true);

  G4Region* defaultRegion = (*(G4RegionStore::GetInstance()))[0];
  RegionInformation* defaultRInfo = new RegionInformation();
  defaultRInfo->SetWorld();
  defaultRInfo->Print();
  defaultRegion->SetUserInformation(defaultRInfo);

  // ------------------------------
  // -- Beampipe
  // ------------------------------
  G4ThreeVector positionBeampipe = G4ThreeVector(0., 0., -10.*cm);
  fSolidBeampipe = new G4Tubs("Beampipe", fBeampipeInnerRadius, fBeampipeOuterRadius, beampipeHalfLength,
			      0.*deg, 360.*deg);
  fLogicBeampipe = new G4LogicalVolume(fSolidBeampipe,fBeampipeMater, "Beampipe", 0, 0, 0);
  fPhysiBeampipe = new G4PVPlacement(0, positionBeampipe, fLogicBeampipe, "Beampipe", fLogicWorld, false, 0, true);
  G4VisAttributes *pBp  = new G4VisAttributes;
  pBp->SetColour(G4Colour(0.7, 0.8, 0.7));
  //  pBp->SetForceLineSegmentsPerCircle(6);
  //  pBp->SetForceAuxEdgeVisible(true);
  pBp->SetLineWidth(0.01);
  pBp->SetForceSolid(false);
  fLogicBeampipe->SetVisAttributes(pBp);

  G4cout << "Beampipe is " << fBeampipeLength/mm << "mm of " << fBeampipeMater->GetName() << G4endl;


  // ------------------------------
  // -- Target
  // ------------------------------
  makeTarget();

  // ------------------------------
  // -- Tracker
  // ------------------------------
  G4ThreeVector positionTracker = G4ThreeVector(0,0,0);

  fSolidTracker = new G4Tubs("Tracker", fTrkInnerRadius, fTrkOuterRadius, trkHalfLength, 0.*deg, 360.*deg);
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


  G4Region* trackerRegion = new G4Region("TrackerRegion");
  RegionInformation* trackerInfo = new RegionInformation();
  trackerInfo->SetTracker();
  trackerRegion->SetUserInformation(trackerInfo);
  fLogicTracker->SetRegion(trackerRegion);
  trackerRegion->AddRootLogicalVolume(fLogicTracker);


  G4double maxStep = 0.5*fChamberWidth;
  fStepLimit = new G4UserLimits(maxStep);
  fLogicTracker->SetUserLimits(fStepLimit);

  // -- shielding
  G4ThreeVector positionShield = G4ThreeVector(0., 0., 70.*cm);
  fSolidShield = new G4Tubs("shielding", 10.0*cm, 150.*cm, 5*cm, 0.*deg, 360.*deg);
  fLogicShield = new G4LogicalVolume(fSolidShield, fAl, "Shield", 0, 0, 0);
  fPhysiShield = new G4PVPlacement(0, positionShield, fLogicShield, "Shield", fLogicWorld, false, 0, true);

  G4VisAttributes *pBs  = new G4VisAttributes;
  pBs->SetColour(G4Colour(0.5, 0.8, 0.5));
  pBs->SetLineWidth(0.01);
  pBs->SetForceSolid(true);

  fLogicShield->SetVisAttributes(pBs);


  // ------------------------------
  // -- Transport tube
  // ------------------------------
  makeSplitTrspTube();
  makeAccel();
  makeEndDetector();


  // ----------------------------------------------------------------------
  // -- Endplates (as a stop for the undecayed Muonia)
  // ----------------------------------------------------------------------
  G4ThreeVector positionEndplate = G4ThreeVector(0., 0., worldHalfLength-0.5*cm);
  fSolidEndplate = new G4Box("Endplate", 200*cm, 200*cm, 0.1*cm);
  fLogicEndplate = new G4LogicalVolume(fSolidEndplate, fTargetMater, "Endplate", 0, 0, 0);
  fPhysiEndplate = new G4PVPlacement(0, positionEndplate, fLogicEndplate, "Endplate", fLogicWorld, false, 0, true);

  fLogicEndplate->SetUserLimits(new G4UserLimits(0.001*mm));

  G4ThreeVector positionEndplateF = G4ThreeVector(0., 0., -(worldHalfLength-0.5*cm));
  fSolidEndplateF = new G4Box("Endplate", 200*cm, 200*cm, 0.1*cm);
  fLogicEndplateF = new G4LogicalVolume(fSolidEndplateF, fTargetMater, "Endplate", 0, 0, 0);
  fPhysiEndplateF = new G4PVPlacement(0, positionEndplateF, fLogicEndplateF, "Endplate", fLogicWorld, false, 0, true);

  fLogicEndplateF->SetUserLimits(new G4UserLimits(0.001*mm));


  G4VisAttributes *pVA1  = new G4VisAttributes;
  pVA1->SetColour(G4Colour(0.9, 0.9, 0.9));
  pVA1->SetForceSolid(true);
  fLogicEndplate->SetVisAttributes(pVA1);
  fLogicEndplateF->SetVisAttributes(pVA1);

  G4cout << "Endplates are " <<  "1mm of " << fTargetMater->GetName()
	 << " with logical name ->" << fLogicEndplate->GetName() << "<-"
	 << G4endl;



  return fPhysiWorld;
}

// ----------------------------------------------------------------------
void DetectorConstruction::SetMagField(G4double fieldValue) {
  G4cout << "===============> DetectorConstruction::SetMagField> SetMagField(" << fieldValue << G4endl;
  fpMagField->SetField(fieldValue);
}


// ----------------------------------------------------------------------
void DetectorConstruction::ConstructSDandField() {

  // -- Sensitive detectors

  // -- FIXME: not sure this is a good thing. Atm I guess that the SD only outputs hits and is not affected by geometry
  if (!G4SDManager::GetSDMpointer()->FindSensitiveDetector("/macs0/TrackerChamberSD", false)) {
    TrackerSD* aTrackerSD = new TrackerSD("/macs0/TrackerChamberSD", "TrackerHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
    SetSensitiveDetector("Chamber_LV", aTrackerSD, true);
  }

  // -- FIXME: not sure this is a good thing. Atm I guess that the SD only outputs hits and is not affected by geometry
  if (!G4SDManager::GetSDMpointer()->FindSensitiveDetector("/macs0/MCPSD", false)) {
    MCPSD* aMCPSD = new MCPSD("/macs0/MCPSD", "MCPHitsCollection");
    G4SDManager::GetSDMpointer()->AddNewDetector(aMCPSD);
    SetSensitiveDetector("lMCP", aMCPSD, true);
  }

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

  G4double orMCP = 2*fTrspOuterRadius;
  G4double lMCP =   0.1*cm;
  auto p0 = new G4Tubs("MCP", 0., orMCP, lMCP, 0.*deg, 360.*deg);
  fLogicMCP = new G4LogicalVolume(p0, fAl, "lMCP");

  G4RotationMatrix* fieldRot = new G4RotationMatrix();
  fieldRot->rotateY(90.*deg);
  fPhysiMCP =  new G4PVPlacement(fieldRot,
				 G4ThreeVector(-(fTrspLength2 + 40.0*cm + 0.5*fTrspOuterRadius),
					       0.,
					       fTrkLength + fTrspLength1 - 0.25*fTrspOuterRadius ),
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

  fTrspOuterRadius = 8.*cm;
  fTrspLength1 = 60.*cm;
  fTrspLength2 = 30.*cm;
  G4double pRtor = 40.*cm;
  fTrspMater  = fAl;
  G4double trkHalfLength  = 0.5*fTrkLength;

  G4ThreeVector positionTrsp = G4ThreeVector(0, 0, 2.*trkHalfLength - 0.25*fTrspLength1);
  G4Tubs  *part1 = new G4Tubs("Trsp1", fTrspOuterRadius-0.1*cm, fTrspOuterRadius, 0.5*fTrspLength1, 0.*deg, 360.*deg);
  G4Tubs  *part3 = new G4Tubs("Trsp3", 2*fTrspOuterRadius-0.1*cm, 2*fTrspOuterRadius, 0.5*fTrspLength2, 0.*deg, 360.*deg);


  // -- trafo for part3
  G4RotationMatrix* rot3 = new G4RotationMatrix;
  rot3->rotateY(M_PI/2.*rad);
  // -- FIXME: the -3*cm are an empirical shift to contain the bending (why? leaking B field?) trajectories inside part3
  G4ThreeVector zTrans3(pRtor + 0.5*fTrspLength2, 0, -0.5*fTrspLength1 - pRtor - 3*cm);
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

// ----------------------------------------------------------------------
void DetectorConstruction::makeTarget() {
  //  G4ThreeVector positionTarget = G4ThreeVector(0., 0., -trkHalfLength+tgtHalfLength);

  G4double zposition = -0.5*fTrkLength + 10.*cm;

  G4ThreeVector positionTarget = G4ThreeVector(0., 0., zposition);
  G4double tgtHalfLength       = 0.5*fTargetLength;    // Half length of the Target

  G4Box* fSolidTarget             = new G4Box("Target", 2*cm, 2*cm, tgtHalfLength);
  G4LogicalVolume* fLogicTarget   = new G4LogicalVolume(fSolidTarget, fTargetMater, "Target", 0, 0, 0);
  G4VPhysicalVolume* fPhysiTarget = new G4PVPlacement(0, positionTarget, fLogicTarget, "Target", fLogicWorld, false, 0, true);

  fLogicTarget->SetUserLimits(new G4UserLimits(0.1*tgtHalfLength));

  G4VisAttributes *pVA  = new G4VisAttributes;
  pVA->SetColour(G4Colour(0.1, 0.8, 0.1));
  pVA->SetForceSolid(true);
  fLogicTarget->SetVisAttributes(pVA);

  G4cout << "Target is " << fTargetLength/mm << "mm of " << fTargetMater->GetName()
	 << " with logical name ->" << fLogicTarget->GetName() << "<-"
	 << G4endl;
}
