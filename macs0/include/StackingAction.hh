#ifndef StackingAction_H
#define StackingAction_H 1

#include "globals.hh"
#include "G4UserStackingAction.hh"

class G4Track;
class G4VHitsCollection;

class StackingAction : public G4UserStackingAction
{
public:
  StackingAction();
  virtual ~StackingAction();

public:
  virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
  virtual void NewStage();
  virtual void PrepareNewEvent();

private:

  G4int fStage;
};

#endif
