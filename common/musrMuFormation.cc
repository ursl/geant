#include "musrMuFormation.hh"

#include "G4ParticleTypes.hh"
#include "RootIO.hh"

using namespace std;

// ----------------------------------------------------------------------
musrMuFormation::musrMuFormation(const G4String& name, G4ProcessType aType)
               : G4VDiscreteProcess(name, aType){

  verboseLevel = 0;
  }

// ----------------------------------------------------------------------
musrMuFormation::~musrMuFormation(){}


// ----------------------------------------------------------------------
G4bool musrMuFormation::IsApplicable(const G4ParticleDefinition& particle) {
  G4cout << "musrMuFormation::IsApplicable for particle " << particle.GetParticleName() << ": "
	 << (&particle == G4MuonPlus::Definition())
	 << G4endl;
  return (&particle == G4MuonPlus::Definition());
}


// ----------------------------------------------------------------------
G4VParticleChange* musrMuFormation::PostStepDoIt(const G4Track& trackData, const G4Step&  aStep) {
  fParticleChange.Initialize(trackData);


  G4Track theNewTrack;
  if (verboseLevel > 0) G4cout << "START verbose = " << verboseLevel
		  << " ID = " << trackData.GetTrackID()
		  << " " << trackData.GetDefinition()->GetParticleName()
		  << " Ekin [keV] = " << trackData.GetKineticEnergy()/CLHEP::keV
		  << " L [mm] = " << trackData.GetTrackLength()/CLHEP::mm
		  << " Edep[keV] = " << aStep.GetTotalEnergyDeposit()/CLHEP::keV
		  << " in " << aStep.GetTrack()->GetVolume()->GetLogicalVolume()->GetName()
		  << G4endl;

  // -- muonium only formed for positive muons entering Target
  G4String  p_name = trackData.GetDefinition()->GetParticleName(); // particle name
  if (p_name != "mu+") {
    if (verboseLevel > 0) G4cout << "musrMuFormation::PostStepDoIt> not a mu, but ->" << p_name << "<-" << G4endl;
    fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;
    return &fParticleChange;
  }
  std::string logVolName = trackData.GetVolume()->GetLogicalVolume()->GetName();
  if (logVolName != "Target") {
    if (verboseLevel > 0) G4cout << "musrMuFormation::PostStepDoIt> not in target, but in ->" << logVolName << "<-" << G4endl;
    fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;
    return &fParticleChange;
  }

  if (verboseLevel > 0) {
    G4cout << "musrMuFormation::PostStepDoIt> mu+ in target, working on Muonium production "
	   << G4endl;
  }

  // --GetDatas(&aStep);
  G4ParticleDefinition *particle;
  G4DynamicParticle *DP;
  G4double yvector[3];
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4double rnd = G4UniformRand();
  G4double E = trackData.GetDynamicParticle()->GetKineticEnergy()/CLHEP::keV;
  fGonin.GetYields(E, 105.658369*1000, yvector); // Energy [keV], muon mass [keV/c2], yield table

  if (verboseLevel > 9) G4cout << "musrMuFormation::GetDatas> rnd = " << rnd
		<< " yvector[0] = " << yvector[0] << " name = " << p_name
		<< G4endl;

  if (rnd < yvector[0]) {
    particle = particleTable->FindParticle(p_name) ;
    if (verboseLevel > 9) G4cout << "musrMuFormation::GetDatas> rnd = " << rnd  << " -> stay with muon"
		  << " yvector[0] = " << yvector[0] << " name = " << p_name
		  << G4endl;
  } else {
    particle = particleTable->FindParticle("Muonium");
    if (verboseLevel > 9) G4cout << "musrMuFormation::GetDatas> rnd = " << rnd << " -> Muonium formed!!!!"
		  << " E = " << E
		  << G4endl
		  << " DP E = " << trackData.GetDynamicParticle()->GetKineticEnergy()
		  << "   vec{|p|} = " << trackData.GetDynamicParticle()->GetMomentumDirection()
		  << G4endl;
  }
  // Set the new dynamic particle DP
  DP = new G4DynamicParticle(particle,
			     trackData.GetDynamicParticle()->GetMomentumDirection(),
			     trackData.GetDynamicParticle()->GetKineticEnergy());
  DP->SetProperTime(  trackData.GetDynamicParticle()->GetProperTime());
  DP->SetPolarization(trackData.GetDynamicParticle()->GetPolarization().x(),
		      trackData.GetDynamicParticle()->GetPolarization().y(),
		      trackData.GetDynamicParticle()->GetPolarization().z());
  DP->SetPreAssignedDecayProperTime(trackData.GetDynamicParticle()->GetPreAssignedDecayProperTime());
  G4Track *aSecondary = new G4Track(DP, trackData.GetGlobalTime(), trackData.GetPosition());

  // -- PrepareSecondary(trackData);
  if (1) {
    if (1) G4cout << "hello" << G4endl;
    RootIO *rio = RootIO::GetInstance();
    TGenVtx *pVtx = rio->getEvent()->addGenVtx();
    pVtx->fNumber = rio->getEvent()->nGenVtx()-1;
    pVtx->fID = -1313;
    pVtx->fvMom.push_back(trackData.GetTrackID());
    pVtx->fvDau.push_back(aSecondary->GetTrackID());
    pVtx->fGlobalTime = trackData.GetGlobalTime();
    pVtx->fLocalTime = trackData.GetLocalTime();
    pVtx->fV.SetX(trackData.GetPosition().x());
    pVtx->fV.SetY(trackData.GetPosition().y());
    pVtx->fV.SetZ(trackData.GetPosition().z());
  }

  fParticleChange.AddSecondary(aSecondary);
  // -- the (original) muon is killed, not the Muonium!
  fParticleChange.ProposeTrackStatus(fStopAndKill) ;

  delete DP;

  return &fParticleChange;
}


// ----------------------------------------------------------------------
G4double musrMuFormation::GetMeanFreePath(const G4Track&t, G4double, G4ForceCondition* condition) {
  *condition = Forced;

  return DBL_MAX;

  std::string logVolName = t.GetVolume()->GetLogicalVolume()->GetName();
  if (verboseLevel > 0) {
    G4cout << "musrMuFormation::GetMeanFreePath> XXXXXXX in " << logVolName
		  << " PDG = " << t.GetDefinition()->GetParticleName()
		  << " ID = " << t.GetTrackID()
	   << G4endl;
  }

  if (logVolName == "Target") {
    return 2*CLHEP::nm;
  }

  return 2*CLHEP::mm;

}
