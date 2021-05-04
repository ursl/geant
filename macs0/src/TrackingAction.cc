#include "TrackingAction.hh"
#include "Trajectory.hh"
#include "TrackInformation.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"

// ----------------------------------------------------------------------
TrackingAction::TrackingAction():G4UserTrackingAction() {
}

// ----------------------------------------------------------------------
void TrackingAction::PreUserTrackingAction(const G4Track* aTrack) {
  TrackInformation* trackInfo = (TrackInformation*)(aTrack->GetUserInformation());
  if (trackInfo->GetTrackingStatus() > 0) {
    // -- these are tracks from primary (1) or secondary (2) particles
    fpTrackingManager->SetStoreTrajectory(true);
    fpTrackingManager->SetTrajectory(new Trajectory(aTrack));
  } else {
    fpTrackingManager->SetStoreTrajectory(false);
  }

}

// ----------------------------------------------------------------------
void TrackingAction::PostUserTrackingAction(const G4Track* aTrack) {
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if (secondaries) {
    TrackInformation* info = (TrackInformation*)(aTrack->GetUserInformation());
    size_t nSeco = secondaries->size();
    if (nSeco > 0) {
      for(size_t i=0; i<nSeco; i++) {
	TrackInformation* infoNew = new TrackInformation(info);
	(*secondaries)[i]->SetUserInformation(infoNew);
      }
    }
  }
}
