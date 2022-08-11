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
// Design of the "new" SMB (as of mid 2022)
// ----------------------------------------------------------------------

using namespace::std;

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
                      fVolume,               //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


  placeSMB();
  
  //  G4AssemblyVolume *solidFibreSMB = makeSMB(fVolume);
  // if (0) {
  //   G4Transform3D TrId;
  //   solidFibreSMB->MakeImprint(fVolume, TrId, 0, 0);  
  // } else {
  // }
  
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
  rotM.rotateY(M_PI*CLHEP::rad);
  rotM.rotateX(M_PI/2*CLHEP::rad);
  rotM.rotateZ(dphi/2);
  G4AssemblyVolume *solidFibreSMB = makeSMB();
  
  for(unsigned int i = 0; i < detector->nribbons; ++i) {
    phi = dphi/2 + i * dphi;
    position =  {-std::sin(phi), std::cos(phi), 0};
    positionPcb = position * (rInSup + rPlate);
    positionPcb.setZ(position.z() - length/2.  - 0.3*CLHEP::cm);
    transform = G4Transform3D(rotM, positionPcb);
    solidFibreSMB->MakeImprint(fVolume, transform);
    rotM.rotateZ(dphi);
  }

  cout << "AAAA> GetAssemblyID() = " << solidFibreSMB->GetAssemblyID() << endl;
  
  vector<G4VPhysicalVolume*>::iterator ipv = solidFibreSMB->GetVolumesIterator();

  while (*ipv) {
    string sname = (*ipv)->GetName(); 
    if (string::npos != sname.find("Pcb")) {
      double phi = (*ipv)->GetTranslation().phi()*57.2957795131; 
      if (phi < 0) phi += 360.;
      int nSmb = 11 - static_cast<int>(phi)/30;
      char ssmb[100];
      sprintf(ssmb, "SMB_US_%d", nSmb);
      cout << "*ipv->GetName() = " << sname << " trsl = "
           << (*ipv)->GetTranslation()
           << " phi = " << phi
           << " -> nSmb = " << nSmb
           << " -> ssmb = " << ssmb
           << endl;
      (*ipv)->SetName(ssmb); 
    }
    ++ipv;
  }
  
}


// ----------------------------------------------------------------------
G4AssemblyVolume* B2DetectorConstruction::makeSMB() {
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
  G4LogicalVolume* fVolumeFibreSMBChip1;
  G4LogicalVolume* fVolumeFibreSMBChip2;
  G4LogicalVolume* fVolumeFibreSMBChip3;
  G4LogicalVolume* fVolumeFibreSMBConnector;
  /*


           +---------------------------------------+
           |                         +--+          |
           |           Chip3         |  |          |
    +------+                         +--+          |
    | C                                            |
    | o                              +--+          |
    | n     C                        |  |          |
    | n      h                       +--+          |
 (2)| e       i                                    | fSMBPcbWidth1
    | c        p                     +--+          |
    | t         1                    |  |          |
    | o                              +--+          |
    | r                                            |
    +------+                         +--+          |
       ^   |         Chip2           |  | (1)      |   (1) fSMBAsicDeltaFront
       |   |                         +--+          |   (2) fSMBPcbWidth2
       |   +---------------------------------------+
       |                fSMBPcbLength1
       |
       +--SMBPcbLength2

  */


  
  double fSMBPcbLength1     = 44.5   * CLHEP::mm;
  double fSMBPcbWidth1      = 26.0   * CLHEP::mm;
  double fSMBPcbThickness   =  1.005 * CLHEP::mm;
  double fSMBPcbLength2     = 11.0   * CLHEP::mm;
  double fSMBPcbWidth2      = 16.0   * CLHEP::mm;

  double fSMBAsicWidth      = 5.0  * CLHEP::mm;
  double fSMBAsicThickness  = 0.3  * CLHEP::mm;

  double fSMBChip1Width     = 5.7  * CLHEP::mm;
  double fSMBChip1Thickness = 0.75 * CLHEP::mm;
  double fSMBChip2Width     = 4.1  * CLHEP::mm;
  double fSMBChip2Thickness = 1.0  * CLHEP::mm;
  double fSMBChip3Width     = 3.1  * CLHEP::mm;
  double fSMBChip3Thickness = 1.0  * CLHEP::mm;

  double fSMBConnectorLength    = 13.1  * CLHEP::mm;
  double fSMBConnectorWidth     =  7.1  * CLHEP::mm;
  double fSMBConnectorThickness =  0.35 * CLHEP::mm;

  double fSMBAsicDeltaFront = 16.3 * CLHEP::mm;    // chip rhs wrt 55.5 (right edge)
  double fSMBAsicDeltaSide  = 1.05 * CLHEP::mm;  
  double fSMBAsicDeltaChip  = 1.30 * CLHEP::mm; 


  double fSMBChip1DeltaCenter = 42.79 * CLHEP::mm; // chip center wrt 55.5 (right edge)
    
  double fSMBChip2DeltaFront = 30.75 * CLHEP::mm;  // chip rhs wrt 55.5 (right edge)
  double fSMBChip2DeltaSide  = 0.25  * CLHEP::mm;

  double fSMBChip3DeltaFront = 28.95 * CLHEP::mm;  // chip rhs wrt 55.5 (right edge) 
  double fSMBChip3DeltaSide  = 22.45 * CLHEP::mm;

  double fSMBConnectorDelta  =  0.355* CLHEP::mm;  // connector separation wrt 0.0 (left edge)

  // -- end header file contents
  
  G4RotationMatrix rotm = G4RotationMatrix();

  // -- the large base
  G4VSolid* solidFibreSMBPcb1   = new G4Box("fibreSMBPcb1",
                                            fSMBPcbWidth1/2., fSMBPcbLength1/2.,
                                            fSMBPcbThickness/2.);
  
  // -- the smaller tail with the connector (interposer?)
  G4VSolid* solidFibreSMBPcb2   = new G4Box("fibreSMBPcb2",
                                            fSMBPcbWidth2/2., fSMBPcbLength2/2.,
                                            fSMBPcbThickness/2.);
  
  // -- the MuTRIG ASIC
  G4VSolid* solidFibreSMBAsic   = new G4Box("fibreSMBAsic",
                                            fSMBAsicWidth/2., fSMBAsicWidth/2.,
                                            fSMBAsicThickness/2.);
  // -- rotated chip next to the interposer
  G4VSolid* solidFibreSMBChip1   = new G4Box("fibreSMBChip1",
                                            fSMBChip1Width/2., fSMBChip1Width/2.,
                                            fSMBChip1Thickness/2.);

  // -- chip with the two round structures *right next* to it
  G4VSolid* solidFibreSMBChip2   = new G4Box("fibreSMBChip2",
                                            fSMBChip2Width/2., fSMBChip2Width/2.,
                                            fSMBChip2Thickness/2.);

  // -- chip with the two round structures at two different edges
  G4VSolid* solidFibreSMBChip3   = new G4Box("fibreSMBChip3",
                                            fSMBChip3Width/2., fSMBChip3Width/2.,
                                            fSMBChip3Thickness/2.);

  // -- Connector at narrow end
  G4VSolid* solidFibreSMBConnector   = new G4Box("fibreSMBConnector",
                                            fSMBConnectorLength/2., fSMBConnectorWidth/2.,
                                            fSMBConnectorThickness/2.);

  
  // -- create PCB (non-rectangular) shape
  G4RotationMatrix* yRot = new G4RotationMatrix(); 
  G4ThreeVector zTrans(0., 0.5*(fSMBPcbLength1+fSMBPcbLength2), 0.);
  G4UnionSolid* solidFibreSMBPcb = new G4UnionSolid("solidFibreSMBPcb",
						    solidFibreSMBPcb1, 
						    solidFibreSMBPcb2,
						    yRot,
						    zTrans);
  
  fVolumeFibreSMBPcb = new G4LogicalVolume(solidFibreSMBPcb,
					   materials.Kapton,
					   "fibreSMBPcb");

  G4VisAttributes *pVA  = new G4VisAttributes;
  pVA->SetColour(G4Colour(0.8, 0.2, 0.4, 0.5));
  pVA->SetForceSolid(true);
  fVolumeFibreSMBPcb->SetVisAttributes(pVA);

  
  // -- create complete board as assembly, first PCB
  G4RotationMatrix Ra;
  G4ThreeVector Ta;
  G4Transform3D Tr;
  G4AssemblyVolume* solidFibreSMB = new G4AssemblyVolume();
  
  //  Ra.rotateZ(M_PI*radian);
  Ra.rotateZ(0.);
  
  Ta.setX(0.);
  Ta.setY(0.);
  Ta.setZ(0.);
  Tr = G4Transform3D(Ra,Ta);
  solidFibreSMB->AddPlacedVolume(fVolumeFibreSMBPcb, Tr);
  
  // -- add 4 asics
  for (unsigned int i = 0; i < 4; ++i) {
    fVolumeFibreSMBAsic[i] = new G4LogicalVolume(solidFibreSMBAsic,
						 materials.Si,
						 "fibreSMBAsic");
    
    G4VisAttributes *pVA2  = new G4VisAttributes;
    pVA2->SetColour(G4Colour(0., 0., 0.));
    pVA2->SetForceSolid(true);
    fVolumeFibreSMBAsic[i]->SetVisAttributes(pVA2);
    
    Ta.setX(0.5*(fSMBPcbWidth1 - fSMBAsicWidth - 2.*fSMBAsicDeltaSide) - i*(fSMBAsicWidth + fSMBAsicDeltaChip));
    Ta.setY(-0.5*(fSMBPcbLength1 - fSMBAsicWidth) + fSMBAsicDeltaFront);
    Ta.setZ(0.5*(fSMBPcbThickness + fSMBAsicThickness));
    
    Tr = G4Transform3D(rotm, Ta);
    solidFibreSMB->AddPlacedVolume(fVolumeFibreSMBAsic[i], Tr);
  }



  
  // -- add chip #2
  fVolumeFibreSMBChip2 = new G4LogicalVolume(solidFibreSMBChip2,
                                             materials.Si,
                                             "fibreSMBChip2");
  
  G4VisAttributes *pVA2  = new G4VisAttributes;
  pVA2->SetColour(G4Colour(0., 0., 0.));
  pVA2->SetForceSolid(true);
  fVolumeFibreSMBChip2->SetVisAttributes(pVA2);
  
  Ta.setX(0.5*(-fSMBPcbWidth1 + fSMBChip2Width + 2.*fSMBChip2DeltaSide));
  Ta.setY(-0.5*(fSMBPcbLength1 - fSMBChip2Width) + fSMBChip2DeltaFront);
  Ta.setZ(0.5*(fSMBPcbThickness + fSMBChip2Thickness));
  
  Tr = G4Transform3D(rotm, Ta);
  solidFibreSMB->AddPlacedVolume(fVolumeFibreSMBChip2, Tr);


  // -- add chip #3
  fVolumeFibreSMBChip3 = new G4LogicalVolume(solidFibreSMBChip3,
                                             materials.Si,
                                             "fibreSMBChip3");
  
  G4VisAttributes *pVA3  = new G4VisAttributes;
  pVA3->SetColour(G4Colour(0., 0., 0.));
  pVA3->SetForceSolid(true);
  fVolumeFibreSMBChip3->SetVisAttributes(pVA3);
  
  Ta.setX(0.5*(-fSMBPcbWidth1 + fSMBChip3Width + 2.*fSMBChip3DeltaSide));
  Ta.setY(-0.5*(fSMBPcbLength1 - fSMBChip3Width) + fSMBChip3DeltaFront);
  Ta.setZ(0.5*(fSMBPcbThickness + fSMBChip3Thickness));
  
  Tr = G4Transform3D(rotm, Ta);
  solidFibreSMB->AddPlacedVolume(fVolumeFibreSMBChip3, Tr);


  // -- add chip #1
  fVolumeFibreSMBChip1 = new G4LogicalVolume(solidFibreSMBChip1,
                                             materials.Si,
                                             "fibreSMBChip1");
  
  G4VisAttributes *pVA1  = new G4VisAttributes;
  pVA1->SetColour(G4Colour(0., 0., 0.));
  pVA1->SetForceSolid(true);
  fVolumeFibreSMBChip1->SetVisAttributes(pVA1);
  
  Ta.setX(0.);
  Ta.setY(-0.5*fSMBPcbLength1 + fSMBChip1DeltaCenter);
  //  Ta.setY(-0.5*(fSMBPcbLength1) + fSMBChip1DeltaFront);
  Ta.setZ(0.5*(fSMBPcbThickness + fSMBChip1Thickness));
  
  rotm.rotateZ(-M_PI/4*CLHEP::rad);
  Tr = G4Transform3D(rotm, Ta);
  solidFibreSMB->AddPlacedVolume(fVolumeFibreSMBChip1, Tr);

  // -- add connector
  fVolumeFibreSMBConnector = new G4LogicalVolume(solidFibreSMBConnector,
                                             materials.Si,
                                             "fibreSMBConnector");
  
  G4VisAttributes *pVA4  = new G4VisAttributes;
  pVA4->SetColour(G4Colour(0.7, 0.7, 0.7));
  pVA4->SetForceSolid(true);
  fVolumeFibreSMBConnector->SetVisAttributes(pVA4);
  
  Ta.setX(0.);
  Ta.setY(0.5*fSMBPcbLength1 + fSMBPcbLength2 - 0.5*fSMBConnectorWidth - fSMBConnectorDelta);
  Ta.setZ(0.5*(fSMBPcbThickness + fSMBConnectorThickness));
  
  // -- rotate back
  rotm.rotateZ(+M_PI/4*CLHEP::rad);
  Tr = G4Transform3D(rotm, Ta);
  solidFibreSMB->AddPlacedVolume(fVolumeFibreSMBConnector, Tr);

  return solidFibreSMB;

  //      int sensorId = (layer << layerOffset) + (ladderId << ladderOffset);
  //      // deliberate offset relative to the grand numbering scheme needed...
  //      ladder->assembly->MakeImprint(mother, trans, sensorId);

}
