#include "musrMuEnergyLossLandau.hh"
#include "musrMuonium.hh"

using namespace std;

// ----------------------------------------------------------------------
musrMuEnergyLossLandau::musrMuEnergyLossLandau(const G4String& name, G4ProcessType aType) :
  G4VDiscreteProcess(name, aType) {
  fRandom = new TRandom3();
  verboseLevel = 0;
  fParticleTable = G4ParticleTable::GetParticleTable();
  fParticleDefinition = fParticleTable->FindParticle("Muonium");
}

// ----------------------------------------------------------------------
musrMuEnergyLossLandau::~musrMuEnergyLossLandau(){}

// double musrMuEnergyLossLandau::landauMPV   = 0.1;
// double musrMuEnergyLossLandau::landauSigma = 0.1;

double musrMuEnergyLossLandau::fLandauMPV   = 0.04;
double musrMuEnergyLossLandau::fLandauSigma = 0.02;


// ----------------------------------------------------------------------
G4bool musrMuEnergyLossLandau::IsApplicable(const G4ParticleDefinition& particle) {
  if (verboseLevel > 0) {
    G4cout << "musrMuEnergyLossLandau::IsApplicable for particle " << particle.GetParticleName() << ": "
	   << (&particle == musrMuonium::Muonium())
	   << G4endl;
  }
  return (&particle == musrMuonium::Muonium());
}

// ----------------------------------------------------------------------
G4VParticleChange* musrMuEnergyLossLandau::PostStepDoIt(const G4Track& trackData, const G4Step&  aStep) {
  fParticleChange.Initialize(trackData);

  // -- check whether Mu hit the endplate. If yes, stop and decay it there.
  G4String p_name = trackData.GetDefinition()->GetParticleName(); // particle name
  if (p_name != "Muonium") {
    return &fParticleChange;
  }

  std::string logVolName = trackData.GetVolume()->GetLogicalVolume()->GetName();
  if (logVolName == "Endplate") {
    if (verboseLevel > 0) {
      G4cout << "musrMuEnergyLossLandau in endplate for  " << p_name
	     << " stopandkill"
	     << G4endl;
    }

    fParticleChange.Initialize(trackData);
    //    fParticleChange.ProposeTrackStatus(fStopAndKill) ;
    fParticleChange.ProposeTrackStatus(fStopButAlive) ;
    return &fParticleChange;
  }

  // -- unwrapped:  if (CheckCondition(aStep)) {
  if (logVolName == "Target") {

    if (verboseLevel > 0) {
      G4cout << "musrMuEnergyLossLandau in target for  " << p_name
	     << " propagate "
	     << G4endl;
    }

    // -- unwrapped:  GetFinalEnergy(&aStep);
    G4double Eloss, Efinal;
    G4double E = trackData.GetDynamicParticle()->GetKineticEnergy()/CLHEP::keV;

    do {
      Eloss = fRandom->Landau(fLandauMPV, fLandauSigma);
      if (Eloss > 0.) break;
    } while (1);

    Efinal = E - Eloss;
    if (Efinal < 0. ) Efinal = 0.;
    if (verboseLevel > 0) G4cout << "musrMuEnergyLossLandau::GetFinalEnergy: "
				 << " E, Eloss, Efinal = " << E << ", " << Eloss << ", "  << Efinal
				 << G4endl;

    // Set the new dynamic particle DP
    DP = new G4DynamicParticle(fParticleDefinition,
			       trackData.GetDynamicParticle()->GetMomentumDirection(),
			       Efinal/1000.);

    DP->SetProperTime(  trackData.GetDynamicParticle()->GetProperTime());
    DP->SetPolarization(trackData.GetDynamicParticle()->GetPolarization().x(),
			trackData.GetDynamicParticle()->GetPolarization().y(),
			trackData.GetDynamicParticle()->GetPolarization().z());
    DP->SetPreAssignedDecayProperTime(trackData.GetDynamicParticle()->GetPreAssignedDecayProperTime());

    //G4Step theStep;
    // NOTE: This is where I should change everything!
    //       - model this according to
    // -- unwrapped:  PrepareSecondary(trackData);
    aSecondary = new G4Track(DP, trackData.GetGlobalTime(), trackData.GetPosition());

    fParticleChange.AddSecondary(aSecondary);
    fParticleChange.ProposeTrackStatus(fStopAndKill) ;
  } else {
    fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;
  }
  return &fParticleChange;
}


// ----------------------------------------------------------------------
G4double musrMuEnergyLossLandau::GetMeanFreePath(const G4Track&t,  G4double, G4ForceCondition* condition1) {
  *condition1 = Forced;
  std::string logVolName = t.GetVolume()->GetLogicalVolume()->GetName();
  if (logVolName == "Target") {
    return 2*CLHEP::nm;
  }
  return 2*CLHEP::mm;
}
