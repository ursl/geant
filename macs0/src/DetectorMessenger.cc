#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "globals.hh"

// ----------------------------------------------------------------------
DetectorMessenger::DetectorMessenger(DetectorConstruction* myDet)
: G4UImessenger(),
  fDetector(myDet),
  fDir(0),
  fDetDir(0),
  fTargMatCmd(0),
  fChamMatCmd(0)
{
  fDir = new G4UIdirectory("/macs0/");
  fDir->SetGuidance("UI commands specific to this example.");

  fDetDir = new G4UIdirectory("/macs0/det/");
  fDetDir->SetGuidance("detector control.");

}

// ----------------------------------------------------------------------
DetectorMessenger::~DetectorMessenger() {
  delete fTargMatCmd;
  delete fChamMatCmd;
  delete fDetDir;
  delete fDir;
}

// ----------------------------------------------------------------------
void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue) {

}
