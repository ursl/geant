#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"

#include "G4GenericMessenger.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "TH1.h"

#include "musrMuonium.hh"


// ----------------------------------------------------------------------
PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* myDC):
  G4VUserPrimaryGeneratorAction(), fVerbose(0), fParticleGun(0), fMyDetector(myDC),
  fSgNpart(0), fBgNpart(1), fBgNpartSigma(0.),
  fSgKinEnergy(1.0), fSgKinEnergySigma(1.1),
  fBgKinEnergy(26.1), fBgKinEnergySigma(1.1) {
  fParticleGun = new G4ParticleGun();

  G4ParticleDefinition* particle  = musrMuonium::MuoniumDefinition();
  fParticleGun->SetParticleDefinition(particle);
  auto particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* backgroundParticle = particleTable->FindParticle("mu+");
  fBackgroundGun = new G4ParticleGun();
  fBackgroundGun->SetParticleDefinition(backgroundParticle);

  DefineCommands();

}

// ----------------------------------------------------------------------
PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fParticleGun;
  delete fBackgroundGun;
}

// ----------------------------------------------------------------------
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  // G4double position = -0.5*(fMyDetector->GetWorldFullLength());
  // fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm, 0.*cm, position));
  // fParticleGun->GeneratePrimaryVertex(anEvent);
  // fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
  // fParticleGun->SetParticleEnergy(1*keV);

  static int first(1);

  if ((fVerbose > 0) || (1 == first)) {
    G4cout << "======================================================================"
	   << G4endl
	   << "==========> Event " << anEvent->GetEventID() << " start."
	   << G4endl
	   << "======================================================================"
	   << G4endl;

    if (fSgNpart > 0) {
      G4cout << "PrimaryGeneratorAction::GeneratePrimaries with particle = "
	     << fParticleGun->GetParticleDefinition()->GetParticleName()
	     << " fSgNpart = " << fSgNpart
	     << G4endl;
    }
  }
  for (int i = 0; i < fSgNpart; ++i) {
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }

  // z0 = -40*cm  is currently at the target position
  G4double z0  = -41.*cm;
  G4double x0  = 0*cm, y0  = 0*cm;
  G4double dx0 = 2*cm, dy0 = 2*cm;
  //  dx0 = 0*cm, dy0 = 0*cm;

  double sig = G4RandGauss::shoot(0., fBgNpartSigma);
  int nBg = round(fBgNpart + sig);

  std::stringstream sstream;
  sstream.setf(std::ios::fixed);
  sstream.precision(6);
  sstream << fBgKinEnergy;

  if ((fVerbose > 0) || (1 == first)) {
    G4cout << "PrimaryGeneratorAction::GeneratePrimaries nbg = " << nBg
	   << Form(" +/- %3.1f", sig)
	   << " Ekin = " << sstream.str()
	   << " MeV at z = " << Form("%5.2f", z0)
	   << Form(" dx,dy =  %3.2f, %3.2f mm", dx0, dy0)
	   << G4endl;
  }

  for (int i = 0; i < nBg; ++i) {
    x0 = dx0*(G4UniformRand()-0.5);
    y0 = dy0*(G4UniformRand()-0.5);

    fBackgroundGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
    fBackgroundGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    fBackgroundGun->SetParticleEnergy(fBgKinEnergy);

    fBackgroundGun->GeneratePrimaryVertex(anEvent);
  }

  // -- turn first false after first call
  if (1 == first) {
    first = 0;
  }
}


// ----------------------------------------------------------------------
void PrimaryGeneratorAction::DefineCommands() {
  fMessenger  = new G4GenericMessenger(this, "/macs0/generator/", "Primary generator control");


  // -- average number of signal particles
  auto& sgNpartCmd = fMessenger->DeclareProperty("sgNpart", fSgNpart,
						 "Average number sg particles.");
  sgNpartCmd.SetParameterName("sgNpart", true);
  sgNpartCmd.SetRange("sgNpart>=0");
  sgNpartCmd.SetDefaultValue("0");


  // -- kinetic energy of signal particles
  auto& sgKinEnergyCmd = fMessenger->DeclarePropertyWithUnit("sgKinEnergy", "MeV", fSgKinEnergy,
							     "Mean kinetic energy of sg particles.");
  sgKinEnergyCmd.SetParameterName("sgEkin", true);
  sgKinEnergyCmd.SetRange("sgEkin>=0.");
  sgKinEnergyCmd.SetDefaultValue("1.");


  // -- sigma of kinetic energy of signal particles
  auto& sgKinEnergySigmaCmd = fMessenger->DeclarePropertyWithUnit("sgKinEnergySigma", "MeV", fSgKinEnergySigma,
								  "Sigma of mean kinetic energy of sg particles.");
  sgKinEnergySigmaCmd.SetParameterName("sgEkinSigma", true);
  sgKinEnergySigmaCmd.SetRange("sgEkinSigma>=0.");
  sgKinEnergySigmaCmd.SetDefaultValue("0.1");



  // -- average number of background particles
  auto& bgNpartCmd = fMessenger->DeclareProperty("bgNpart", fBgNpart,
						 "Average number bg particles.");
  bgNpartCmd.SetParameterName("bgNpart", true);
  bgNpartCmd.SetRange("bgNpart>=0");
  bgNpartCmd.SetDefaultValue("2");

  // -- sigma of average number of background particles
  auto& bgNpartSigmaCmd = fMessenger->DeclareProperty("bgNpartSigma", fBgNpartSigma,
						 "Sigma of average number bg particles.");
  bgNpartSigmaCmd.SetParameterName("bgNpartSigma", true);
  bgNpartSigmaCmd.SetRange("bgNpartSigma>=0.");
  bgNpartSigmaCmd.SetDefaultValue("0.2");


  // -- kinetic energy of background particles
  auto& bgKinEnergyCmd = fMessenger->DeclarePropertyWithUnit("bgKinEnergy", "MeV", fBgKinEnergy,
							     "Mean kinetic energy of bg particles.");
  bgKinEnergyCmd.SetParameterName("bgEkin", true);
  bgKinEnergyCmd.SetRange("bgEkin>=0.");
  bgKinEnergyCmd.SetDefaultValue("26.");


  // -- sigma of kinetic energy of background particles
  auto& bgKinEnergySigmaCmd = fMessenger->DeclarePropertyWithUnit("bgKinEnergySigma", "MeV", fBgKinEnergySigma,
								  "Sigma of mean kinetic energy of bg particles.");
  bgKinEnergySigmaCmd.SetParameterName("bgEkinSigma", true);
  bgKinEnergySigmaCmd.SetRange("bgEkinSigma>=0.");
  bgKinEnergySigmaCmd.SetDefaultValue("1.");

}
