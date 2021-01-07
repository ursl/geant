#include "RunAction.hh"
#include "G4Run.hh"
#include "RootIO.hh"

// ----------------------------------------------------------------------
RunAction::RunAction() : G4UserRunAction() {}

// ----------------------------------------------------------------------
RunAction::~RunAction() {}

// ----------------------------------------------------------------------
void RunAction::BeginOfRunAction(const G4Run* aRun) {
  G4cout << "==========> Run " << aRun->GetRunID() << " start." << G4endl;
  RootIO *rio = RootIO::GetInstance();
  rio->getEvent()->fRunNumber = aRun->GetRunID();
}

// ----------------------------------------------------------------------
void RunAction::EndOfRunAction(const G4Run*) { }
