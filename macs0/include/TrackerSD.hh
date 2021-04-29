#ifndef TrackerSD_h
#define TrackerSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TrackerHit.hh"

class G4Step;
class G4HCofThisEvent;

// ----------------------------------------------------------------------
class TrackerSD : public G4VSensitiveDetector {
  public:
  TrackerSD(const G4String& name, const G4String& hitsCollectionName, int verbose = 0);
  ~TrackerSD();

  virtual void Initialize(G4HCofThisEvent*);
  virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent*);

private:
  TrackerHitsCollection* fHitsCollection;
  int fVerbose;
};

#endif
