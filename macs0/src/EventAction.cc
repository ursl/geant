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
  //  G4cout << "==========> Event " << evt->GetEventID() << " start." << G4endl;
}

// ----------------------------------------------------------------------
void EventAction::EndOfEventAction(const G4Event* evt) {
  G4int event_id = evt->GetEventID();
  RootIO *rio = RootIO::GetInstance();
  rio->getEvent()->fEventNumber = event_id;

  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  bool DBX(true);
  double px, py, pz, m, ekin, etot, vx, vy, vz;
  int pdgid, pid, tid;

  if (DBX) G4cout << "----------------------------------------------------------------------" << G4endl;
  for (G4int i = 0; i < n_trajectories; ++i) {
    G4Trajectory* trj=(G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
    ekin = trj->GetInitialKineticEnergy();

    tid   = trj->GetTrackID();
    pid   = trj->GetParentID();
    m     = trj->GetParticleDefinition()->GetPDGMass();
    pdgid = trj->GetPDGEncoding();

    px    = trj->GetInitialMomentum().x();
    py    = trj->GetInitialMomentum().y();
    pz    = trj->GetInitialMomentum().z();
    vx    = trj->GetPoint(0)->GetPosition().x();
    vy    = trj->GetPoint(0)->GetPosition().y();
    vz    = trj->GetPoint(0)->GetPosition().z();

    ekin  = trj->GetInitialKineticEnergy();
    etot  = TMath::Sqrt(px*px + py*py + pz*pz + m*m);
    if (DBX) G4cout <<
	       Form("%3d ID= %+5d trkId= %2d parentID= %2d (E, p)= %+7.2f/%+7.3f/%+7.3f/%+7.3f ekin= %7.4f v=(%+4.3f, %+4.3f, %+13.8f)",
		    i, pdgid, tid, pid, etot,
		    px, py, pz, ekin,
		    vx, vy, vz
		    )
		    << G4endl;

    TGenCand* pGen =  rio->getEvent()->addGenCand();
    pGen->fID = pdgid;
    pGen->fNumber = trj->GetTrackID();
    pGen->fMom1   = trj->GetParentID();
    pGen->fMom2   = -1;
    pGen->fDau1   = 9999;
    pGen->fDau2   = -9999;
    pGen->fStatus = 0;
    pGen->fQ      = trj->GetCharge();
    pGen->fP.SetXYZT(px, py, pz, etot);
    pGen->fV.SetXYZ(trj->GetPoint(0)->GetPosition().x(),
		    trj->GetPoint(0)->GetPosition().y(),
		    trj->GetPoint(0)->GetPosition().z());
  }

  // -- now determine the daughter index pointers
  TGenCand *pMom(0), *pDau(0);
  for (int imom = 0; imom < rio->getEvent()->nGenCands(); ++imom) {
    pMom = rio->getEvent()->getGenCand(imom);
    for (int idau = imom; idau < rio->getEvent()->nGenCands(); ++idau) {
      pDau = rio->getEvent()->getGenCand(idau);
      if (pMom->fNumber == pDau->fMom1) {
	if (pDau->fMom1 < pMom->fDau1 ) {
	  pMom->fDau1 = pDau->fNumber;
	}
	if (pDau->fMom1 > pMom->fDau2 ) {
	  pMom->fDau2 = pDau->fNumber;
	}
      }
    }
  }

  rio->fillTree();

  rio->getEvent()->Clear();
}
