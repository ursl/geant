#include "MCPHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<MCPHit> MCPHitAllocator;
//G4ThreadLocal G4Allocator<MCPHit>* MCPHitAllocator = 0;

// ----------------------------------------------------------------------
MCPHit::MCPHit(): G4VHit(), fTrackID(0), fChamberNb(0), fEdep(0), fPos(0,0,0) {
}

// ----------------------------------------------------------------------
MCPHit::~MCPHit() {
}


// ----------------------------------------------------------------------
MCPHit::MCPHit(const MCPHit& right) : G4VHit() {
  fTrackID   = right.fTrackID;
  fChamberNb = right.fChamberNb;
  fEdep      = right.fEdep;
  fPos       = right.fPos;
}

// ----------------------------------------------------------------------
const MCPHit& MCPHit::operator=(const MCPHit& right) {
  fTrackID   = right.fTrackID;
  fChamberNb = right.fChamberNb;
  fEdep      = right.fEdep;
  fPos       = right.fPos;
  return *this;
}

// ----------------------------------------------------------------------
G4int MCPHit::operator==(const MCPHit& right) const {
  return (this==&right) ? 1 : 0;
}

// ----------------------------------------------------------------------
void MCPHit::Draw() {
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager) {
    G4cout << "=============> MCPHit::Draw()> Draw hit at " << fPos << G4endl;
    G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(0.8,0.8,0.6);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  } else {
    G4cout << "=============> MCPHit::Draw()> NO G4VVisManager found!!" << G4endl;
  }

}

// ----------------------------------------------------------------------
void MCPHit::Print() {
  G4cout << "  MCP trackID: " << fTrackID << "  chNb: " << fChamberNb
         << "  Edep[MeV]: " << fEdep
         << "  pos[mm]: " << fPos
	 << "  Etrk[MeV]: " << fEtrk
	 << G4endl;
}
