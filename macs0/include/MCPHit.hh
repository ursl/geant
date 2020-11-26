#ifndef MCPHit_h
#define MCPHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

// ----------------------------------------------------------------------
class MCPHit : public G4VHit {
public:

  MCPHit();
  ~MCPHit();
  MCPHit(const MCPHit&);
  const MCPHit& operator=(const MCPHit&);
  G4int operator==(const MCPHit&) const;

  inline void* operator new(size_t);
  inline void  operator delete(void*);

  virtual void Draw();
  virtual void Print();

public:

  void SetTrackID  (G4int track)      { fTrackID = track; };
  void SetChamberNb(G4int chamb)      { fChamberNb = chamb; };
  void SetEdep     (G4double de)      { fEdep = de; };
  void SetEtrk     (G4double e)       { fEtrk = e; };
  void SetGblTime  (G4double t)       { fGblTime = t; };
  void SetPos      (G4ThreeVector xyz){ fPos = xyz; };

  G4int GetTrackID()    { return fTrackID; };
  G4int GetChamberNb()  { return fChamberNb; };
  G4double GetEdep()    { return fEdep; };
  G4double GetEtrk()    { return fEtrk; };
  G4double GetGblTime() { return fGblTime; };
  G4ThreeVector GetPos(){ return fPos; };

private:

  G4int         fTrackID;
  G4int         fChamberNb;
  G4double      fEdep;
  G4double      fEtrk;
  G4double      fGblTime;
  G4ThreeVector fPos;
};


typedef G4THitsCollection<MCPHit> MCPHitsCollection;

extern G4Allocator<MCPHit> MCPHitAllocator;
//extern G4ThreadLocal G4Allocator<MCPHit>* MCPHitAllocator;


// ----------------------------------------------------------------------
inline void* MCPHit::operator new(size_t) {
  void *aHit;
  aHit = (void *) MCPHitAllocator.MallocSingle();
  return aHit;
}


// ----------------------------------------------------------------------
inline void MCPHit::operator delete(void *aHit) {
  MCPHitAllocator.FreeSingle((MCPHit*) aHit);
}


#endif
