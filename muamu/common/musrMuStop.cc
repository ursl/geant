#include "musrMuStop.hh"

using namespace std;

musrMuStop::musrMuStop(const G4String& name, G4ProcessType  aType) : G4VDiscreteProcess(name, aType) {}

musrMuStop:: ~musrMuStop(){}

// ----------------------------------------------------------------------
G4VParticleChange* musrMuStop::PostStepDoIt(const G4Track& trackData, const G4Step& aStep) {
  fParticleChange.Initialize(trackData);

  //! Tao - Get time information */
  itime = trackData.GetProperTime();
  gtime = trackData.GetGlobalTime();
  ftime = trackData.GetDynamicParticle()->GetPreAssignedDecayProperTime();

  deltatime = ftime - itime;
  fParticleChange.ProposeGlobalTime(deltatime + itime -gtime);

  /*! - Set position, momentum, energy and time of the particle change. */
  fParticleChange.ProposePosition(trackData.GetPosition());
  fParticleChange.ProposeMomentumDirection(trackData.GetMomentumDirection());
  fParticleChange.ProposeEnergy(trackData.GetKineticEnergy());
  fParticleChange.ProposeGlobalTime(gtime);
  fParticleChange.ProposeProperTime(itime);
  fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;

  //  if (trackData.GetKineticEnergy() <

  if (CheckCondition(aStep)) {
    fParticleChange.ProposePosition(trackData.GetStep()->GetPreStepPoint()->GetPosition());
    fParticleChange.ProposeTrackStatus(fStopButAlive) ;
  }

  /*! - Return the changed particle object. */
  return &fParticleChange;
}


// ----------------------------------------------------------------------
// -- Muonium will be stopped as soon as it enters a material different than vacuum or target.
G4bool musrMuStop::CheckCondition(const G4Step& aStep) {
  G4bool condition = false;
  p_name = aStep.GetTrack()->GetDefinition()->GetParticleName(); // particle name
  if (p_name == "Muonium"
      && aStep.GetTrack()->GetVolume()->GetLogicalVolume()->GetName() != "Target"
      && aStep.GetTrack()->GetVolume()->GetLogicalVolume()->GetMaterial()->GetName() != "vac"
      ) {
    condition=true;
    //    condition = false; //ul: do not stop
    }
  return condition;
}

// ----------------------------------------------------------------------
G4double musrMuStop::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;
}

// ----------------------------------------------------------------------
void musrMuStop::PrepareSecondary(const G4Track& track) {
  aSecondary = new G4Track(DP,track.GetDynamicParticle()->GetPreAssignedDecayProperTime(),track.GetPosition());
}
