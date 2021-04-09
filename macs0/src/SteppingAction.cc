#include "SteppingAction.hh"
#include "RegionInformation.hh"
#include "TrackInformation.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Region.hh"
#include "G4SteppingManager.hh"

// ----------------------------------------------------------------------
SteppingAction::SteppingAction(): G4UserSteppingAction(){

}

// ----------------------------------------------------------------------
SteppingAction::~SteppingAction() {
}

// ----------------------------------------------------------------------
void SteppingAction::UserSteppingAction(const G4Step * theStep) {
  // Suspend a track if it is entering into the calorimeter

  // check if it is alive
  G4Track * theTrack = theStep->GetTrack();
  if (theTrack->GetTrackStatus() != fAlive) {
    return;
  }

  // get region information
  G4StepPoint * thePrePoint       = theStep->GetPreStepPoint();
  G4LogicalVolume * thePreLV      = thePrePoint->GetPhysicalVolume()->GetLogicalVolume();
  RegionInformation* thePreRInfo  = (RegionInformation*)(thePreLV->GetRegion()->GetUserInformation());
  G4StepPoint * thePostPoint      = theStep->GetPostStepPoint();
  G4LogicalVolume * thePostLV     = thePostPoint->GetPhysicalVolume()->GetLogicalVolume();
  RegionInformation* thePostRInfo = (RegionInformation*)(thePostLV->GetRegion()->GetUserInformation());

  // check if it is entering to the calorimeter volume
  if (!(thePreRInfo->IsCalorimeter()) && (thePostRInfo->IsCalorimeter()))  {
    // if the track had already been suspended at the previous step, let it go.
    TrackInformation* trackInfo  = static_cast<TrackInformation*>(theTrack->GetUserInformation());
    if (trackInfo->GetSuspendedStepID()>-1) {
      if (fpSteppingManager->GetverboseLevel() > 0) {
        G4cout<<"++++ This track had already been suspended at step #"
	      <<trackInfo->GetSuspendedStepID()<<". Tracking resumed."
	      <<G4endl;
      }
    } else {
      trackInfo->SetSuspendedStepID(theTrack->GetCurrentStepNumber());
      theTrack->SetTrackStatus(fSuspend);
      if (fpSteppingManager->GetverboseLevel() > 0) {
	G4cout<<"++++ This track now being suspended" << G4endl;
      }
    }
  }
}
