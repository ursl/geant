#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

class G4GenericMessenger;
class DetectorConstruction;
class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGeneratorAction(DetectorConstruction*);
  ~PrimaryGeneratorAction();

public:
  virtual void GeneratePrimaries(G4Event*);
  void setVerbose(int v) {fVerbose = v;}

private:
  int      fVerbose;
  G4int    fSgNpart;
  G4int    fBgNpart;
  G4double fBgNpartSigma;
  G4double fSgKinEnergy, fSgKinEnergySigma;
  G4double fBgKinEnergy, fBgKinEnergySigma;

  G4double fSgGunZposition;
  G4double fBgGunZposition;

  G4ParticleGun* fParticleGun;
  G4ParticleGun* fBackgroundGun;
  DetectorConstruction* fMyDetector;

  void DefineCommands();
  G4GenericMessenger* fMessenger;
};

#endif
