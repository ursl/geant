#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

#include "RootIO.hh"

// ----------------------------------------------------------------------
EventAction::EventAction(): G4UserEventAction() { }

// ----------------------------------------------------------------------
EventAction::~EventAction() {
}

// ----------------------------------------------------------------------
void EventAction::BeginOfEventAction(const G4Event*evt) {
  G4cout << "==========> Event " << evt->GetEventID() << " start." << G4endl;
}

// ----------------------------------------------------------------------
void EventAction::EndOfEventAction(const G4Event* evt) {
  G4int event_id = evt->GetEventID();
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  RootIO *rio = RootIO::GetInstance();
  std::vector<int> genlist;
  G4cout << "----------------------------------------------------------------------" << G4endl;
  for (G4int i = 0; i < n_trajectories; ++i) {
    G4Trajectory* trj=(G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
    if (0 == trj->GetParentID()) {
      genlist.push_back(i);
    }
    if (1 == trj->GetParentID()) {
      genlist.push_back(i);
    }
    G4cout << " Evt/trj =  " << event_id << "/" << i << " " << trj->GetParticleName()
	   << " ID/Mum " << trj->GetTrackID() << "/" << trj->GetParentID()
	   << " charge: " << trj->GetCharge()
	   << " energy: " << trj->GetInitialKineticEnergy()
	   <<  G4endl;
  }

  G4VHitsCollection* hc = evt->GetHCofThisEvent()->GetHC(0);
  G4cout << " Hits stored in this event:  "  << hc->GetSize() << G4endl;

  G4cout << "----------------------------------------------------------------------" << G4endl;
  G4cout << "genlist size() = " << genlist.size() << G4endl;

  std::string sid;
  double px, py, pz, e;
  int pdgid, pid, tid;
  for (unsigned int i = 0; i < genlist.size(); ++i) {
    G4Trajectory* trj=(G4Trajectory*)((*(evt->GetTrajectoryContainer()))[genlist[i]]);
    px = trj->GetInitialMomentum().x();
    py = trj->GetInitialMomentum().y();
    pz = trj->GetInitialMomentum().z();
    e  = trj->GetInitialKineticEnergy();
    pdgid = trj->GetPDGEncoding();
    tid = trj->GetTrackID();
    pid = trj->GetParentID();
    if (pdgid == 11) sid = "e-";
    if (pdgid == -11) sid = "e+";
    if (pdgid == 12) sid = "nu_e";
    if (pdgid == -12) sid = "anti-nu_e";
    if (pdgid == 13) sid = "mu-";
    if (pdgid == -13) sid = "mu+";
    if (pdgid == 14) sid = "nu_mu";
    if (pdgid == -14) sid = "anti-nu_mu";
    G4cout << "(p,Ekin) = " << px << "/" << py << "/" << pz << "/" << e << " ID = " << pdgid
	   << " (" << sid << ", m = " << trj->GetParticleDefinition()->GetPDGMass() << ")"
	   << " trkID = " << tid << " parentID = " << pid
	   << G4endl;
    TGenCand* pGen =  rio->getEvent()->addGenCand();
    pGen->fID = pdgid;
    pGen->fNumber = genlist[i];
    pGen->fP.SetXYZT(px, py, px, e);
  }

  rio->fillTree();
  rio->getEvent()->Clear();
}
