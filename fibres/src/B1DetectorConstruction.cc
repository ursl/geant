#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"


#include <G4UnionSolid.hh>

// ----------------------------------------------------------------------
B1DetectorConstruction::B1DetectorConstruction(): G4VUserDetectorConstruction(), fScoringVolume(0) { }


// ----------------------------------------------------------------------
B1DetectorConstruction::~B1DetectorConstruction() { }


// ----------------------------------------------------------------------
G4VPhysicalVolume* B1DetectorConstruction::Construct() {  
  G4NistManager* nist = G4NistManager::Instance();
  
  // -- Envelope parameters
  G4double env_sizeXY = 20*cm, env_sizeZ = 30*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_Air");
   
  // -- Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  // -- World
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld = new G4Box("World",              
				0.5*world_sizeXY,
				0.5*world_sizeXY,
				0.5*world_sizeZ);     
      
  G4LogicalVolume* volume = new G4LogicalVolume(solidWorld,
						world_mat,        
						"World"); 
  
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      volume,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  // -- begin placeholder material
  struct material {
    G4Material* Si;  
    G4Material* SiO2;
    G4Material *Kapton, *PVC;
  } materials;
  materials.Si     = nist->FindOrBuildMaterial("G4_SI");
  materials.SiO2   = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  materials.Kapton = nist->FindOrBuildMaterial("G4_KAPTON");
  materials.PVC    = nist->FindOrBuildMaterial("G4_PVC");
  
  // -- begin header file contents
  G4LogicalVolume* fVolumeFibreFEE;
  G4LogicalVolume* fVolumeFibreFEEPcb;
  G4LogicalVolume* fVolumeFibreFEEAsic;

  
  double fFEEPcbLength1    = 43.5 * CLHEP::mm;
  double fFEEPcbWidth1     = 25.6 * CLHEP::mm;
  double fFEEPcbThickness  = 1.07 * CLHEP::mm;
  double fFEEPcbLength2    = 11.1 * CLHEP::mm;
  double fFEEPcbWidth2     = 15.7 * CLHEP::mm;
  double fFEEAsicLength    = 5.0  * CLHEP::mm;
  double fFEEAsicWidth     = 5.0  * CLHEP::mm;
  double fFEEAsicThickness = 1.0  * CLHEP::mm; // FIXME??
  // -- end header file contents
  
  
  G4VSolid* solidFibreFEEPcb   = new G4Box("fibreFEEPcb",
					   fFEEPcbWidth1/2., fFEEPcbLength1/2.,
					   fFEEPcbThickness/2.);
 
  G4VSolid* solidFibreFEEAsic   = new G4Box("fibreFEEAsic",
					    fFEEAsicWidth/2., fFEEAsicLength/2.,
					    fFEEAsicThickness/2.);
  
  
  G4RotationMatrix* yRot = new G4RotationMatrix(); 
  G4ThreeVector zTrans(0, 10, 10);
  G4UnionSolid* solidFibreFEE = new G4UnionSolid("solidFibreFEEPcb+solidFibreFEEAsic",
						 solidFibreFEEPcb,
						 solidFibreFEEAsic,
						 yRot,
						 zTrans);
  
  fVolumeFibreFEEPcb = new G4LogicalVolume(solidFibreFEEPcb,
					   materials.PVC,
					   "fibreFEEPcb");
  
  fVolumeFibreFEE = new G4LogicalVolume(solidFibreFEE,
					materials.Si,
					"fibreFEE");
  

  double lengthPlate = 20*mm;
  double length      = 10*mm;
  
  new G4PVPlacement(0, {0, 0, -lengthPlate/2. + length/2.},
		    fVolumeFibreFEEPcb,
		    "fibreFEEPcb",
		    volume,
		    false,
		    0);

  return physWorld;
}

