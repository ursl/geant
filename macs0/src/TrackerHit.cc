#include "TrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<TrackerHit> TrackerHitAllocator;

// ----------------------------------------------------------------------
TrackerHit::TrackerHit(): G4VHit(), fTrackID(0), fChamberNb(0), fEdep(0), fPos(0,0,0) {
}

// ----------------------------------------------------------------------
TrackerHit::~TrackerHit() {
}


// ----------------------------------------------------------------------
TrackerHit::TrackerHit(const TrackerHit& right) : G4VHit() {
  fTrackID   = right.fTrackID;
  fChamberNb = right.fChamberNb;
  fEdep      = right.fEdep;
  fPos       = right.fPos;
}

// ----------------------------------------------------------------------
const TrackerHit& TrackerHit::operator=(const TrackerHit& right) {
  fTrackID   = right.fTrackID;
  fChamberNb = right.fChamberNb;
  fEdep      = right.fEdep;
  fPos       = right.fPos;
  return *this;
}

// ----------------------------------------------------------------------
G4int TrackerHit::operator==(const TrackerHit& right) const {
  return (this==&right) ? 1 : 0;
}

// ----------------------------------------------------------------------
void TrackerHit::Draw() {
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager) {
    G4Circle circle(fPos);
    circle.SetScreenSize(2.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

// ----------------------------------------------------------------------
void TrackerHit::Print() {
  G4cout << "  trackID: " << fTrackID << "  chamberNb: " << fChamberNb
         << "  energy deposit[MeV]: " << fEdep
         << "  position[mm]: " << fPos << G4endl;
}
