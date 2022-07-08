#include "B3DetectorConstruction.hh"

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

#include "G4VisAttributes.hh"

#include <G4AssemblyVolume.hh>
#include <G4UnionSolid.hh>
#include <G4MultiUnion.hh>
#include <G4Tubs.hh>
#include <G4Polycone.hh>
#include <G4SubtractionSolid.hh>

// ----------------------------------------------------------------------
B3DetectorConstruction::B3DetectorConstruction(): G4VUserDetectorConstruction(), fScoringVolume(0) { }


// ----------------------------------------------------------------------
B3DetectorConstruction::~B3DetectorConstruction() { }


// ----------------------------------------------------------------------
G4VPhysicalVolume* B3DetectorConstruction::Construct() {  
  G4NistManager* nist = G4NistManager::Instance();

  // -- Envelope parameters
  G4double env_sizeXY = 200*cm, env_sizeZ = 300*cm;
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
      
  fVolume = new G4LogicalVolume(solidWorld,
						world_mat,        
						"World"); 
  
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      fVolume,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  G4Polycone* solidFibreColdFlange = makeColdFlange(fVolume);
  G4LogicalVolume* logicFibreColdFlange =                         
    new G4LogicalVolume(solidFibreColdFlange,//its solid
                        env_mat,             //its material
                        "ColdFlange");       //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicFibreColdFlange,    //its logical volume
                    "ColdFlange",            //its name
                    fVolume,                 //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  //  solidFibreColdFlange->MakeImprint(fVolume, transform);
  
  return physWorld;
}



// ----------------------------------------------------------------------
G4Polycone* B3DetectorConstruction::makeColdFlange(G4LogicalVolume *volume) {
  G4NistManager* nist = G4NistManager::Instance();

  // -- begin placeholder material
  struct material {
    G4Material* Si;  
    G4Material* SiO2;
    G4Material *Kapton, *PVC;
  } materials;
  materials.Si = nist->FindOrBuildMaterial("G4_SI");
  materials.SiO2 = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  materials.Kapton = nist->FindOrBuildMaterial("G4_KAPTON");
  materials.PVC = nist->FindOrBuildMaterial("G4_PVC");
  // -- end placeholder material
  
  // -- begin header file contents
  double phiStart = M_PI/4.*CLHEP::rad,
    phiTotal = 3/2*M_PI*CLHEP::rad;

  int numZPlanes = 9;
  double rInner[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
  double rOuter[] = { 0, 10, 10, 5 , 5, 10 , 10 , 2, 2};
  double zPlane[] = { 5, 7, 9, 11, 25, 27, 29, 31, 35 };

  G4Polycone* solidFibreColdFlange = new G4Polycone("coldFlange",
                                                          phiStart,
                                                          phiTotal,
                                                          numZPlanes,
                                                          zPlane,
                                                          rInner,
                                                          rOuter);
  
  return solidFibreColdFlange;
}
