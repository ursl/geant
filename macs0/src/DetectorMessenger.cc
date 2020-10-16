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
  fN02Dir(0),
  fDetDir(0),
  fTargMatCmd(0),
  fChamMatCmd(0),
  fFieldCmd(0)
{
  fN02Dir = new G4UIdirectory("/muamu/");
  fN02Dir->SetGuidance("UI commands specific to this example.");

  fDetDir = new G4UIdirectory("/muamu/det/");
  fDetDir->SetGuidance("detector control.");

  fFieldCmd = new G4UIcmdWithADoubleAndUnit("/muamu/det/setField",this);
  fFieldCmd->SetGuidance("Define magnetic field.");
  fFieldCmd->SetGuidance("Magnetic field will be in X direction.");
  fFieldCmd->SetParameterName("Bx",false);
  fFieldCmd->SetUnitCategory("Magnetic flux density");
  fFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

// ----------------------------------------------------------------------
DetectorMessenger::~DetectorMessenger() {
  delete fTargMatCmd;
  delete fChamMatCmd;
  delete fFieldCmd;
  delete fDetDir;
  delete fN02Dir;
}

// ----------------------------------------------------------------------
void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue) {
  if (command == fFieldCmd) { fDetector->SetMagField(fFieldCmd->GetNewDoubleValue(newValue));}
}
