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

#include "G4VisAttributes.hh"

#include <G4AssemblyVolume.hh>
#include <G4UnionSolid.hh>
#include <G4MultiUnion.hh>
#include <G4Tubs.hh>
#include <G4SubtractionSolid.hh>

// ----------------------------------------------------------------------
B1DetectorConstruction::B1DetectorConstruction(): G4VUserDetectorConstruction(), fScoringVolume(0) { }


// ----------------------------------------------------------------------
B1DetectorConstruction::~B1DetectorConstruction() { }


// ----------------------------------------------------------------------
G4VPhysicalVolume* B1DetectorConstruction::Construct() {  
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


  G4AssemblyVolume *solidFibreFEE = makeFEE(fVolume);
  if (1) {
    G4Transform3D TrId;
    solidFibreFEE->MakeImprint(fVolume, TrId, 0, 0);  
  } else {
    placeFEE();
  }
  
  return physWorld;
}


// ----------------------------------------------------------------------
void B1DetectorConstruction::placeFEE() {
  struct bla {
    int nribbons;
  };
  struct bla *detector = new struct bla;
  detector->nribbons = 12;
  double rInSup =   39 * CLHEP::mm;
  double rPlate =   15 * CLHEP::mm;
  double length = (20.66 * mm) * 18/2. + (0.04 * mm) * (18 - 1) + 2; 
    
  G4RotationMatrix rotM = G4RotationMatrix();
  G4ThreeVector position = {};
  G4ThreeVector positionPcb = {};
  G4Transform3D transform;
  
  double phi;
  const double dphi = 2 * M_PI / detector->nribbons;
  // -- y = 0 is between two ribbons. Therefore start with offset. Mirror this here as well.
  rotM.rotateX(-M_PI/2*CLHEP::rad);
  rotM.rotateZ(0);
  G4AssemblyVolume *solidFibreFEE = makeFEE(fVolume);
  
  for(unsigned int i = 0; i < detector->nribbons; ++i) {
    //	phi = dphi/2 + i * dphi;
    phi = i * dphi;
    position =  {-std::sin(phi), std::cos(phi), 0};
    positionPcb = position * (rInSup + rPlate);
    positionPcb.setZ(position.z() - length/2.  - 0.3*CLHEP::cm);
    transform = G4Transform3D(rotM, positionPcb);
    solidFibreFEE->MakeImprint(fVolume, transform);
    rotM.rotateZ(dphi);
  }
}


// ----------------------------------------------------------------------
G4AssemblyVolume* B1DetectorConstruction::makeFEE(G4LogicalVolume *volume) {
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
  G4LogicalVolume* fVolumeFibreFEE;
  G4LogicalVolume* fVolumeFibreFEEPcb;
  G4LogicalVolume* fVolumeFibreFEEAsic[4];

  
  double fFEEPcbLength1     = 43.5 * CLHEP::mm;
  double fFEEPcbWidth1      = 25.6 * CLHEP::mm;
  double fFEEPcbThickness   = 1.07 * CLHEP::mm;
  double fFEEPcbLength2     = 11.1 * CLHEP::mm;
  double fFEEPcbWidth2      = 15.7 * CLHEP::mm;
  double fFEEAsicLength     = 5.0  * CLHEP::mm;
  double fFEEAsicWidth      = 5.0  * CLHEP::mm;
  double fFEEAsicThickness  = 1.0  * CLHEP::mm; // FIXME??
  double fFEEAsicDeltaFront = 7.9  * CLHEP::mm; 
  double fFEEAsicDeltaSide  = 2*1.236* CLHEP::mm; 
  double fFEEAsicDeltaChip  = 1.05 * CLHEP::mm; 
  // -- end header file contents
  
  
  G4VSolid* solidFibreFEEPcb1   = new G4Box("fibreFEEPcb1",
					    fFEEPcbWidth1/2., fFEEPcbLength1/2.,
					    fFEEPcbThickness/2.);

  G4MultiUnion* fourHoles = new G4MultiUnion("FourHoles");
  G4Tubs *hole = new G4Tubs("hole", 0, 1*mm, fFEEPcbThickness, 0, 2.*M_PI*radian);
  G4RotationMatrix rotm = G4RotationMatrix();
  G4Transform3D tr1 = G4Transform3D(rotm, G4ThreeVector(-11, -19, 0.));
  fourHoles->AddNode(*hole, tr1);
  G4Transform3D tr2 = G4Transform3D(rotm, G4ThreeVector(-11, 19, 0.));
  fourHoles->AddNode(*hole, tr2);
  G4Transform3D tr3 = G4Transform3D(rotm, G4ThreeVector(+11, -19, 0.));
  fourHoles->AddNode(*hole, tr3);
  G4Transform3D tr4 = G4Transform3D(rotm, G4ThreeVector(+11, 19, 0.));
  fourHoles->AddNode(*hole, tr4);
  
  G4SubtractionSolid *subtraction = new G4SubtractionSolid("fibreFEEPcb11", solidFibreFEEPcb1, fourHoles);

  
  G4VSolid* solidFibreFEEPcb2   = new G4Box("fibreFEEPcb2",
					    fFEEPcbWidth2/2., fFEEPcbLength2/2.,
					    fFEEPcbThickness/2.);
  
  G4VSolid* solidFibreFEEAsic   = new G4Box("fibreFEEAsic",
					    fFEEAsicWidth/2., fFEEAsicLength/2.,
					    fFEEAsicThickness/2.);
  
  
  // -- create PCB (non-rectangular) shape
  G4RotationMatrix* yRot = new G4RotationMatrix(); 
  G4ThreeVector zTrans(0., 0.5*(fFEEPcbLength1+fFEEPcbLength2), 0.);
  G4UnionSolid* solidFibreFEEPcb = new G4UnionSolid("solidFibreFEEPcb",
						    subtraction, 
						    solidFibreFEEPcb2,
						    yRot,
						    zTrans);
  
  fVolumeFibreFEEPcb = new G4LogicalVolume(solidFibreFEEPcb,
					   materials.Kapton,
					   "fibreFEEPcb");

  G4VisAttributes *pVA1  = new G4VisAttributes;
  pVA1->SetColour(G4Colour(0.8, 0.2, 0.4, 0.5));
  pVA1->SetForceSolid(true);
  fVolumeFibreFEEPcb->SetVisAttributes(pVA1);

  
  // -- create complete board as assembly, first PCB
  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;
  G4AssemblyVolume* solidFibreFEE = new G4AssemblyVolume();
  
  Ra.rotateZ(M_PI*radian);
  
  Ta.setX(0.);
  Ta.setY(0.);
  Ta.setZ(0.);
  Tr = G4Transform3D(Ra,Ta);
  solidFibreFEE->AddPlacedVolume(fVolumeFibreFEEPcb, Tr);
  
  // -- add asics
  for (unsigned int i = 0; i < 4; ++i) {
    fVolumeFibreFEEAsic[i] = new G4LogicalVolume(solidFibreFEEAsic,
						 materials.Si,
						 "fibreFEEAsic");
    
    G4VisAttributes *pVA2  = new G4VisAttributes;
    pVA2->SetColour(G4Colour(0., 0., 0.));
    pVA2->SetForceSolid(true);
    fVolumeFibreFEEAsic[i]->SetVisAttributes(pVA2);
    
    //    Ta.setX(-0.5*fFEEPcbWidth1 + 0.5*(fFEEAsicDeltaSide + fFEEAsicWidth) + i*(fFEEAsicWidth + fFEEAsicDeltaChip)); 
    Ta.setX(0.5*(-fFEEPcbWidth1 + fFEEAsicDeltaSide + fFEEAsicWidth) + i*(fFEEAsicWidth + fFEEAsicDeltaChip));
    Ta.setY(0.5*fFEEPcbLength1 - 0.5*fFEEAsicWidth - fFEEAsicDeltaFront);
    Ta.setZ(0.5*(fFEEPcbThickness + fFEEAsicThickness));
    
    Tr = G4Transform3D(rotm, Ta);
    solidFibreFEE->AddPlacedVolume(fVolumeFibreFEEAsic[i], Tr);
  }

  return solidFibreFEE;
}
