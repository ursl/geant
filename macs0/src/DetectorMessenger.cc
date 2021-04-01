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
  fTargLenCmd(0),
  fModMatCmd(0),
  fModLenCmd(0)
{
  fDir = new G4UIdirectory("/macs0/");
  fDir->SetGuidance("UI commands specific to this example.");

  fDetDir = new G4UIdirectory("/macs0/det/");
  fDetDir->SetGuidance("detector control.");


  fTargMatCmd = new G4UIcmdWithAString("/macs0/det/setTargetMaterial", this);
  fTargMatCmd->SetGuidance("Select Material of the Target.");
  fTargMatCmd->SetParameterName("choice", false);
  fTargMatCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fTargLenCmd = new G4UIcmdWithADoubleAndUnit("/macs0/det/setTargetLength", this);
  fTargLenCmd->SetGuidance("Set length (thickness) of target");
  fTargLenCmd->SetParameterName("length", false);
  fTargLenCmd->SetRange("length>0.");
  fTargLenCmd->SetUnitCategory("Length");
  fTargLenCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  fTargLenCmd->SetToBeBroadcasted(false);

}

// ----------------------------------------------------------------------
DetectorMessenger::~DetectorMessenger() {
  delete fTargMatCmd;
  delete fDetDir;
  delete fDir;
}

// ----------------------------------------------------------------------
void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {
  if (command == fTargMatCmd) {
    G4cout << "==DetectorMessenger::SetNewValue> calling SetTargetMaterial"
	   << G4endl;
    fDetector->SetTargetMaterial(newValue);
  }

  if (command == fTargLenCmd) {
    G4cout << "==DetectorMessenger::SetNewValue> calling SetTargetLength"
	   << G4endl;
    fDetector->SetTargetLength(fTargLenCmd->GetNewDoubleValue(newValue));
  }


}
