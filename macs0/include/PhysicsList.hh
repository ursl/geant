#ifndef PhysicsList_h
#define PhysicsList_h 1

#include <list>

#include "G4VUserPhysicsList.hh"

class G4VPhysicsConstructor;
class G4BertiniElectroNuclearBuilder;

class PhysicsList: public G4VUserPhysicsList {
public:
  PhysicsList();
  virtual ~PhysicsList();
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  virtual void SetCuts();

private:
  std::list<G4VPhysicsConstructor*> fPhysicsConstructors;
  G4BertiniElectroNuclearBuilder* bertiniElectroNuclearBuilder = nullptr;

  G4double fEGammaCut, fMuonPolarization;

};

#endif
