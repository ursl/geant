#ifndef TrackerHit_h
#define TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

// ----------------------------------------------------------------------
class TrackerHit : public G4VHit {
public:

  TrackerHit();
  ~TrackerHit();
  TrackerHit(const TrackerHit&);
  const TrackerHit& operator=(const TrackerHit&);
  G4int operator==(const TrackerHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  virtual void Draw();
  virtual void Print();

public:

  void SetTrackID  (G4int track)      { fTrackID = track; };
  void SetChamberNb(G4int chamb)      { fChamberNb = chamb; };
  void SetEdep     (G4double de)      { fEdep = de; };
  void SetPos      (G4ThreeVector xyz){ fPos = xyz; };

  G4int GetTrackID()    { return fTrackID; };
  G4int GetChamberNb()  { return fChamberNb; };
  G4double GetEdep()    { return fEdep; };
  G4ThreeVector GetPos(){ return fPos; };

private:

  G4int         fTrackID;
  G4int         fChamberNb;
  G4double      fEdep;
  G4ThreeVector fPos;
};


typedef G4THitsCollection<TrackerHit> TrackerHitsCollection;

extern G4Allocator<TrackerHit> TrackerHitAllocator;


// ----------------------------------------------------------------------
inline void* TrackerHit::operator new(size_t) {
  void *aHit;
  aHit = (void *) TrackerHitAllocator.MallocSingle();
  return aHit;
}


// ----------------------------------------------------------------------
inline void TrackerHit::operator delete(void *aHit) {
  TrackerHitAllocator.FreeSingle((TrackerHit*) aHit);
}


#endif
