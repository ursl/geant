#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;


class RunAction : public G4UserRunAction {
public:
  RunAction();
  ~RunAction();

public:
  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);
  void setVerbose(int v) {fVerbose = v;}

private:
  int fVerbose;

};

#endif
