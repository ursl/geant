#include "StackingAction.hh"
#include "TrackInformation.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

// ----------------------------------------------------------------------
StackingAction::StackingAction():G4UserStackingAction(), fStage(0){
}

// ----------------------------------------------------------------------
StackingAction::~StackingAction() {
}

// ----------------------------------------------------------------------
G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track * aTrack) {
  G4ClassificationOfNewTrack classification = fUrgent;
  if (0 == fStage) {
    TrackInformation* trackInfo;
    // Track reached calorimeter
    if (fSuspend == aTrack->GetTrackStatus())   {
      trackInfo = (TrackInformation*)(aTrack->GetUserInformation());
      trackInfo->SetTrackingStatus(0);
      trackInfo->SetSourceTrackInformation(aTrack);
      if (0) G4cout << "Track " << aTrack->GetTrackID()
		    << " (parentID " << aTrack->GetParentID()
		    << ") has reached calorimeter and has been suspended at "
		    << aTrack->GetPosition() << G4endl;
      classification = fWaiting;
    } else if (0 == aTrack->GetParentID())  {
      if (0) G4cout << "Track " << aTrack->GetTrackID()  << " from primary particle " << G4endl;
      // Primary particle
      trackInfo = new TrackInformation(aTrack);
      trackInfo->SetTrackingStatus(1);
      G4Track* theTrack = (G4Track*)aTrack;
      theTrack->SetUserInformation(trackInfo);
    } else {
      if (0) G4cout << "Track " << aTrack->GetTrackID()  << " from other particle " << G4endl;
      // secondary particle
      trackInfo = new TrackInformation(aTrack);
      trackInfo->SetTrackingStatus(2);
      G4Track* theTrack = (G4Track*)aTrack;
      theTrack->SetUserInformation(trackInfo);
    }
  }
  return classification;
}


// ----------------------------------------------------------------------
void StackingAction::NewStage() {
  // G4cout << "+++++++++++ Stage " << fStage << G4endl;
  if (0 == fStage)  {
    // display trajetory information in the tracking region
    G4cout << G4endl;
    G4cout << "Tracks in tracking region have been processed. -- Stage 0 over."
           << G4endl;
    G4cout << G4endl;
  } else  {
    G4cout << "StackingAction::NewStage> fStage != 0" << G4endl;
  }

  if (stackManager->GetNUrgentTrack()) {
    // Transfer all tracks in Urgent stack to Waiting stack, since all tracks
    // in Waiting stack have already been transfered to Urgent stack before
    // invokation of this method.
    stackManager->TransferStackedTracks(fUrgent,fWaiting);

    // Then, transfer only one track to Urgent stack.
    stackManager->TransferOneStackedTrack(fWaiting,fUrgent);

    fStage++;
  }
}

// ----------------------------------------------------------------------
void StackingAction::PrepareNewEvent() {
  fStage = 0;
}
