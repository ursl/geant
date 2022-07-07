#include "B2DetectorConstruction.hh"

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
B2DetectorConstruction::B2DetectorConstruction(): G4VUserDetectorConstruction(), fScoringVolume(0) { }


// ----------------------------------------------------------------------
B2DetectorConstruction::~B2DetectorConstruction() { }


// ----------------------------------------------------------------------
G4VPhysicalVolume* B2DetectorConstruction::Construct() {  
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


  G4AssemblyVolume *solidFibreSMB = makeSMB(fVolume);
  if (0) {
    G4Transform3D TrId;
    solidFibreSMB->MakeImprint(fVolume, TrId, 0, 0);  
  } else {
    placeSMB();
  }
  
  return physWorld;
}


// ----------------------------------------------------------------------
void B2DetectorConstruction::placeSMB() {
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
  rotM.rotateZ(dphi/2);
  G4AssemblyVolume *solidFibreSMB = makeSMB(fVolume);
  
  for(unsigned int i = 0; i < detector->nribbons; ++i) {
    phi = dphi/2 + i * dphi;
    position =  {-std::sin(phi), std::cos(phi), 0};
    positionPcb = position * (rInSup + rPlate);
    positionPcb.setZ(position.z() - length/2.  - 0.3*CLHEP::cm);
    transform = G4Transform3D(rotM, positionPcb);
    solidFibreSMB->MakeImprint(fVolume, transform);
    rotM.rotateZ(dphi);
  }
}


// ----------------------------------------------------------------------
G4AssemblyVolume* B2DetectorConstruction::makeSMB(G4LogicalVolume *volume) {
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
  G4LogicalVolume* fVolumeFibreSMB;
  G4LogicalVolume* fVolumeFibreSMBPcb;
  G4LogicalVolume* fVolumeFibreSMBAsic[4];

  
  double fSMBPcbLength1     = 43.5 * CLHEP::mm;
  double fSMBPcbWidth1      = 25.6 * CLHEP::mm;
  double fSMBPcbThickness   = 1.07 * CLHEP::mm;
  double fSMBPcbLength2     = 11.1 * CLHEP::mm;
  double fSMBPcbWidth2      = 15.7 * CLHEP::mm;
  double fSMBAsicLength     = 5.0  * CLHEP::mm;
  double fSMBAsicWidth      = 5.0  * CLHEP::mm;
  double fSMBAsicThickness  = 1.0  * CLHEP::mm; // FIXME??
  double fSMBAsicDeltaFront = 7.9  * CLHEP::mm; 
  double fSMBAsicDeltaSide  = 2*1.236* CLHEP::mm;  // need factor 2 to better center ASIC row
  double fSMBAsicDeltaChip  = 1.05 * CLHEP::mm; 
  // -- end header file contents
  
  
  G4VSolid* solidFibreSMBPcb1   = new G4Box("fibreSMBPcb1",
					    fSMBPcbWidth1/2., fSMBPcbLength1/2.,
					    fSMBPcbThickness/2.);

  G4MultiUnion* fourHoles = new G4MultiUnion("FourHoles");
  G4Tubs *hole = new G4Tubs("hole", 0, 1*mm, fSMBPcbThickness, 0, 2.*M_PI*radian);
  G4RotationMatrix rotm = G4RotationMatrix();
  G4Transform3D tr1 = G4Transform3D(rotm, G4ThreeVector(-11, -19, 0.));
  fourHoles->AddNode(*hole, tr1);
  G4Transform3D tr2 = G4Transform3D(rotm, G4ThreeVector(-11, 19, 0.));
  fourHoles->AddNode(*hole, tr2);
  G4Transform3D tr3 = G4Transform3D(rotm, G4ThreeVector(+11, -19, 0.));
  fourHoles->AddNode(*hole, tr3);
  G4Transform3D tr4 = G4Transform3D(rotm, G4ThreeVector(+11, 19, 0.));
  fourHoles->AddNode(*hole, tr4);
  
  G4SubtractionSolid *subtraction = new G4SubtractionSolid("fibreSMBPcb11", solidFibreSMBPcb1, fourHoles);

  
  G4VSolid* solidFibreSMBPcb2   = new G4Box("fibreSMBPcb2",
					    fSMBPcbWidth2/2., fSMBPcbLength2/2.,
					    fSMBPcbThickness/2.);
  
  G4VSolid* solidFibreSMBAsic   = new G4Box("fibreSMBAsic",
					    fSMBAsicWidth/2., fSMBAsicLength/2.,
					    fSMBAsicThickness/2.);
  
  
  // -- create PCB (non-rectangular) shape
  G4RotationMatrix* yRot = new G4RotationMatrix(); 
  G4ThreeVector zTrans(0., 0.5*(fSMBPcbLength1+fSMBPcbLength2), 0.);
  G4UnionSolid* solidFibreSMBPcb = new G4UnionSolid("solidFibreSMBPcb",
						    subtraction, 
						    solidFibreSMBPcb2,
						    yRot,
						    zTrans);
  
  fVolumeFibreSMBPcb = new G4LogicalVolume(solidFibreSMBPcb,
					   materials.Kapton,
					   "fibreSMBPcb");

  G4VisAttributes *pVA1  = new G4VisAttributes;
  pVA1->SetColour(G4Colour(0.8, 0.2, 0.4, 0.5));
  pVA1->SetForceSolid(true);
  fVolumeFibreSMBPcb->SetVisAttributes(pVA1);

  
  // -- create complete board as assembly, first PCB
  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;
  G4AssemblyVolume* solidFibreSMB = new G4AssemblyVolume();
  
  Ra.rotateZ(M_PI*radian);
  
  Ta.setX(0.);
  Ta.setY(0.);
  Ta.setZ(0.);
  Tr = G4Transform3D(Ra,Ta);
  solidFibreSMB->AddPlacedVolume(fVolumeFibreSMBPcb, Tr);
  
  // -- add asics
  for (unsigned int i = 0; i < 4; ++i) {
    fVolumeFibreSMBAsic[i] = new G4LogicalVolume(solidFibreSMBAsic,
						 materials.Si,
						 "fibreSMBAsic");
    
    G4VisAttributes *pVA2  = new G4VisAttributes;
    pVA2->SetColour(G4Colour(0., 0., 0.));
    pVA2->SetForceSolid(true);
    fVolumeFibreSMBAsic[i]->SetVisAttributes(pVA2);
    
    //    Ta.setX(-0.5*fSMBPcbWidth1 + 0.5*(fSMBAsicDeltaSide + fSMBAsicWidth) + i*(fSMBAsicWidth + fSMBAsicDeltaChip)); 
    Ta.setX(0.5*(-fSMBPcbWidth1 + fSMBAsicDeltaSide + fSMBAsicWidth) + i*(fSMBAsicWidth + fSMBAsicDeltaChip));
    Ta.setY(0.5*fSMBPcbLength1 - 0.5*fSMBAsicWidth - fSMBAsicDeltaFront);
    Ta.setZ(0.5*(fSMBPcbThickness + fSMBAsicThickness));
    
    Tr = G4Transform3D(rotm, Ta);
    solidFibreSMB->AddPlacedVolume(fVolumeFibreSMBAsic[i], Tr);
  }

  return solidFibreSMB;
}
