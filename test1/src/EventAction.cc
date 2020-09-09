#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

// ----------------------------------------------------------------------
EventAction::EventAction(): G4UserEventAction() { }

// ----------------------------------------------------------------------
EventAction::~EventAction() { }

// ----------------------------------------------------------------------
void EventAction::BeginOfEventAction(const G4Event*) { }

// ----------------------------------------------------------------------
void EventAction::EndOfEventAction(const G4Event* evt) {
  G4int event_id = evt->GetEventID();

  // get number of stored trajectories
  //
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  if (event_id < 100 || event_id%100 == 0) {
    G4cout << ">>> Event " << evt->GetEventID() << G4endl;
    G4cout << "    " << n_trajectories
           << " trajectories stored in this event." << G4endl;
  }
}
