#include "TrackerHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<TrackerHit> TrackerHitAllocator;
//G4ThreadLocal G4Allocator<TrackerHit>* TrackerHitAllocator = 0;

// ----------------------------------------------------------------------
TrackerHit::TrackerHit(): G4VHit(), fTrackID(0), fChamberNb(0), fEdep(0), fGblTime(0), fPos(0,0,0) {
}

// ----------------------------------------------------------------------
TrackerHit::~TrackerHit() {
}


// ----------------------------------------------------------------------
TrackerHit::TrackerHit(const TrackerHit& right) : G4VHit() {
  fTrackID   = right.fTrackID;
  fChamberNb = right.fChamberNb;
  fEdep      = right.fEdep;
  fEtrk      = right.fEtrk;
  fPos       = right.fPos;
  fGblTime   = right.fGblTime;
}

// ----------------------------------------------------------------------
const TrackerHit& TrackerHit::operator=(const TrackerHit& right) {
  fTrackID   = right.fTrackID;
  fChamberNb = right.fChamberNb;
  fEdep      = right.fEdep;
  fPos       = right.fPos;
  fEtrk      = right.fEtrk;
  fGblTime   = right.fGblTime;
  return *this;
}

// ----------------------------------------------------------------------
G4int TrackerHit::operator==(const TrackerHit& right) const {
  return (this==&right) ? 1 : 0;
}

// ----------------------------------------------------------------------
void TrackerHit::Draw() {
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (pVVisManager) {
    G4cout << "=============> TrackerHit::Draw()> Draw hit at " << fPos << G4endl;
    G4Circle circle(fPos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(0.4,0.8,0.6);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  } else {
    G4cout << "=============> TrackerHit::Draw()> NO G4VVisManager found!!" << G4endl;
  }

}

// ----------------------------------------------------------------------
void TrackerHit::Print() {
  G4cout << "  TRK trackID: " << fTrackID << "  chamberNb: " << fChamberNb
         << "  energy deposit[MeV]: " << fEdep
         << "  position[mm]: " << fPos
	 << " Ekin: " << fEtrk
	 << " GblTime: " << fGblTime
	 << G4endl;
}
