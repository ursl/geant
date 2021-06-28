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
  fSgNpart(0), fBgNpart(0), fBgNpartSigma(0.),
  fSgKinEnergy(1.0), fSgKinEnergySigma(1.1), fSgAlpha(0.),
  fBgKinEnergy(26.1), fBgKinEnergySigma(1.1) {
  fParticleGun = new G4ParticleGun();
  fSgGunZposition = -39.*cm;;
  fBgGunZposition = -41.*cm;;

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

  G4double psiMin = 0*deg;       //psi in [0, 2*pi]
  G4double psiMax = 360*deg;
  G4double alphaMax = fSgAlpha;
  G4double cosAlphaMin = std::cos(0.);
  G4double cosAlphaMax = std::cos(alphaMax);

  G4double z0  = fSgGunZposition;
  G4double x0  = 0*cm, y0  = 0*cm;
  G4double dx0 = 2*cm, dy0 = 2*cm;

  if ((fVerbose > 0) || (1 == first)) {
    std::stringstream sstream;
    sstream.setf(std::ios::fixed);
    sstream.precision(6);
    sstream << fSgKinEnergy;
    G4cout << "======================================================================"
	   << G4endl
	   << "==========> PrimaryGeneratorAction::GeneratePrimaries " << anEvent->GetEventID() << " start."
	   << G4endl
	   << "======================================================================"
	   << G4endl;

    if (fSgNpart > 0) {
      G4cout << "PrimaryGeneratorAction::GeneratePrimaries with particle = "
	     << fParticleGun->GetParticleDefinition()->GetParticleName()
	     << " fSgNpart = " << fSgNpart
	     << " Ekin = " << sstream.str()
	     << " MeV at z = " << Form("%5.2f", z0)
	     << Form(" dx,dy =  %3.2f, %3.2f mm", dx0, dy0)
	     << " alpha (max) = " << Form("%3.1f", fSgAlpha/deg) << " deg"
	     << G4endl;
    }
  }

  fParticleGun->SetParticleEnergy(fSgKinEnergy);
  for (int i = 0; i < fSgNpart; ++i) {
    x0 = dx0*(G4UniformRand() - 0.5);
    y0 = dy0*(G4UniformRand() - 0.5);

    fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));

    if (fSgAlpha > 0.) {
      // examples//extended/eventgenerator/particleGun/src/PrimaryGeneratorAction0.cc
      G4double cosAlpha = cosAlphaMin-G4UniformRand()*(cosAlphaMin - cosAlphaMax);
      G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);
      G4double psi = psiMin + G4UniformRand()*(psiMax - psiMin);

      G4double ux = sinAlpha*std::cos(psi),
	uy = sinAlpha*std::sin(psi),
	uz = cosAlpha;
      G4ThreeVector dir(ux,uy,uz);

      fParticleGun->SetParticleMomentumDirection(dir);
      G4cout << "Mu direction: " << dir
	     << " alpha = " << dir.angle(G4ThreeVector(0., 0., 1.))/deg
	     << G4endl;
    } else {
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    }


    fParticleGun->GeneratePrimaryVertex(anEvent);
  }


  // -- now background
  z0  = fBgGunZposition;
  x0  = 0*cm;
  y0  = 0*cm;
  dx0 = 2*cm;
  dy0 = 2*cm;

  double sig = G4RandGauss::shoot(0., fBgNpartSigma);
  int nBg = round(fBgNpart + sig);

  if ((fVerbose > 0) || (1 == first)) {
    std::stringstream sstream;
    sstream.setf(std::ios::fixed);
    sstream.precision(6);
    sstream << fBgKinEnergy;
    G4cout << "PrimaryGeneratorAction::GeneratePrimaries nbg = " << nBg
	   << Form(" +/- %3.1f", sig)
	   << " Ekin = " << sstream.str()
	   << " MeV at z = " << Form("%5.2f", z0)
	   << Form(" dx,dy =  %3.2f, %3.2f mm", dx0, dy0)
	   << G4endl;
  }

  fBackgroundGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  fBackgroundGun->SetParticleEnergy(fBgKinEnergy);
  for (int i = 0; i < nBg; ++i) {
    x0 = dx0*(G4UniformRand()-0.5);
    y0 = dy0*(G4UniformRand()-0.5);

    fBackgroundGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
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

  // -- z position of sg gun
  auto& sgGunZpositionCmd = fMessenger->DeclarePropertyWithUnit("sgGunZposition", "cm", fSgGunZposition,
								  "z position of sg particle gun.");
  sgGunZpositionCmd.SetParameterName("sgGunZposition", true);
  sgGunZpositionCmd.SetDefaultValue("-39.");

  // -- opening angle for signal particles
  auto& sgAlphaCmd = fMessenger->DeclarePropertyWithUnit("sgAlpha", "deg", fSgAlpha,
								  "Maximum opening angle for signal particles.");
  sgAlphaCmd.SetParameterName("sgAlpha", true);
  sgAlphaCmd.SetRange("sgAlpha>=0.");
  sgAlphaCmd.SetDefaultValue("0.");



  // -- average number of background particles
  auto& bgNpartCmd = fMessenger->DeclareProperty("bgNpart", fBgNpart,
						 "Average number bg particles.");
  bgNpartCmd.SetParameterName("bgNpart", true);
  bgNpartCmd.SetRange("bgNpart>=0");
  bgNpartCmd.SetDefaultValue("0");

  // -- sigma of average number of background particles
  auto& bgNpartSigmaCmd = fMessenger->DeclareProperty("bgNpartSigma", fBgNpartSigma,
						 "Sigma of average number bg particles.");
  bgNpartSigmaCmd.SetParameterName("bgNpartSigma", true);
  bgNpartSigmaCmd.SetRange("bgNpartSigma>=0.");
  bgNpartSigmaCmd.SetDefaultValue("0.0");


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

  // -- z position of bg gun
  auto& bgGunZpositionCmd = fMessenger->DeclarePropertyWithUnit("bgGunZposition", "cm", fBgGunZposition,
								  "z position of bg particle gun.");
  bgGunZpositionCmd.SetParameterName("bgGunZposition", true);
  bgGunZpositionCmd.SetDefaultValue("-41.");

}
