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
	   << " Ekin: " << trj->GetInitialKineticEnergy()
	   <<  G4endl;
  }

  rio->getEvent()->fEventNumber = event_id;

  G4VHitsCollection* hc = evt->GetHCofThisEvent()->GetHC(0);

  std::string sid;
  double px, py, pz, m, ekin, etot;
  int pdgid, pid, tid, nchg(0);
  for (unsigned int i = 0; i < genlist.size(); ++i) {
    G4Trajectory* trj=(G4Trajectory*)((*(evt->GetTrajectoryContainer()))[genlist[i]]);
    tid   = trj->GetTrackID();
    pid   = trj->GetParentID();
    m     = trj->GetParticleDefinition()->GetPDGMass();
    pdgid = trj->GetPDGEncoding();
    if (pdgid == 11) sid = "e-";
    if (pdgid == -11) sid = "e+";
    if (pdgid == 12) sid = "nu_e";
    if (pdgid == -12) sid = "anti-nu_e";
    if (pdgid == 13) sid = "mu-";
    if (pdgid == -13) sid = "mu+";
    if (pdgid == 14) sid = "nu_mu";
    if (pdgid == -14) sid = "anti-nu_mu";
    if (pdgid == 22) sid = "gamma";
    if (11 == TMath::Abs(pdgid) || 13 == TMath::Abs(pdgid)) ++nchg;

    px    = trj->GetInitialMomentum().x();
    py    = trj->GetInitialMomentum().y();
    pz    = trj->GetInitialMomentum().z();
    ekin  = trj->GetInitialKineticEnergy();
    etot  = TMath::Sqrt(px*px + py*py + pz*pz + m*m);
    G4cout << Form("(E, p) = %5.4f/%5.4f/%5.4f/%5.4f ID = %4d(%s) m = %5.4f trkId = %2d parentID = %2d",
		   etot, px, py, pz, pdgid, sid.c_str(), trj->GetParticleDefinition()->GetPDGMass(),
		   tid, pid)
	   << G4endl;
    TGenCand* pGen =  rio->getEvent()->addGenCand();
    pGen->fID = pdgid;
    pGen->fNumber = genlist[i];
    pGen->fP.SetXYZT(px, py, pz, etot);
  }

  if (nchg > 1) rio->fillTree();

  rio->getEvent()->Clear();
}
