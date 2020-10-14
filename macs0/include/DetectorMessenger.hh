#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class DetectorMessenger: public G4UImessenger {
public:
  DetectorMessenger(DetectorConstruction*);
  ~DetectorMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);

private:
  DetectorConstruction* fDetector;

  G4UIdirectory*             fN02Dir;
  G4UIdirectory*             fDetDir;
  G4UIcmdWithAString*        fTargMatCmd;
  G4UIcmdWithAString*        fChamMatCmd;
  G4UIcmdWithADoubleAndUnit* fFieldCmd;
};

#endif
