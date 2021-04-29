#ifndef MCPSD_h
#define MCPSD_h 1

#include "G4VSensitiveDetector.hh"
#include "MCPHit.hh"

class G4Step;
class G4HCofThisEvent;

// ----------------------------------------------------------------------
class MCPSD : public G4VSensitiveDetector {
  public:
  MCPSD(const G4String& name, const G4String& hitsCollectionName, int verbose = 0);
  ~MCPSD();

  virtual void Initialize(G4HCofThisEvent*);
  virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  virtual void EndOfEvent(G4HCofThisEvent*);

private:
  MCPHitsCollection* fHitsCollection;
  int fVerbose;
};

#endif
