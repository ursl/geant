#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

#include "Trajectory.hh"

#include "G4ios.hh"

#include "RootIO.hh"
#include <TTimeStamp.h>

// ----------------------------------------------------------------------
EventAction::EventAction(): G4UserEventAction(), fVerbose(0) { }

// ----------------------------------------------------------------------
EventAction::~EventAction() {
}

// ----------------------------------------------------------------------
void EventAction::BeginOfEventAction(const G4Event*evt) {
  TTimeStamp ts;
  int every(1);
  if (0 == fVerbose) {
    if (fNevt <= 100) {
      every = 10;
    } else if (fNevt <= 1000) {
      every = 100;
    } else if (fNevt <= 10000) {
      every = 500;
    } else {
      every = 1000;
    }
  }
  if (0 == evt->GetEventID()%every) {
    G4cout << "==========> Event " << Form("%6d", evt->GetEventID()) << " start, time now: "
	   << ts.AsString("lc")
	   << G4endl;
  }
}

// ----------------------------------------------------------------------
void EventAction::EndOfEventAction(const G4Event* evt) {
  G4int event_id = evt->GetEventID();
  RootIO *rio = RootIO::GetInstance();
  rio->getEvent()->fEventNumber = event_id;

  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  bool DBX(fVerbose>0);
  double px, py, pz, m, ekin, etot, vx, vy, vz;
  int pdgid, parentIdx, tIdx;

  if (DBX) G4cout << "----------------------------------------------------------------------" << G4endl;
  for (G4int i = 0; i < n_trajectories; ++i) {
    Trajectory* trj = (Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
    ekin = trj->GetKineticEnergy();

    tIdx   = trj->GetTrackID();
    parentIdx   = trj->GetParentID();
    m     = trj->GetMass();
    pdgid = trj->GetPDGEncoding();

    px    = trj->GetInitialMomentum().x();
    py    = trj->GetInitialMomentum().y();
    pz    = trj->GetInitialMomentum().z();
    vx    = trj->GetPoint(0)->GetPosition().x();
    vy    = trj->GetPoint(0)->GetPosition().y();
    vz    = trj->GetPoint(0)->GetPosition().z();

    etot  = TMath::Sqrt(px*px + py*py + pz*pz + m*m);

    TGenCand* pGen =  rio->getEvent()->addGenCand();
    pGen->fID = pdgid;
    pGen->fNumber = trj->GetTrackID();
    pGen->fMom1   = trj->GetParentID();
    pGen->fMom2   = -1;
    pGen->fDau1   = 9999;
    pGen->fDau2   = -9999;
    pGen->fStatus = 0;
    pGen->fQ      = trj->GetCharge();
    pGen->fMass   = m;
    pGen->fLocalTime  = trj->GetLocalTime();
    pGen->fGlobalTime = trj->GetGlobalTime();
    pGen->fP.SetXYZT(px, py, pz, etot);
    pGen->fV.SetXYZ(trj->GetPoint(0)->GetPosition().x(),
		    trj->GetPoint(0)->GetPosition().y(),
		    trj->GetPoint(0)->GetPosition().z());

    if (DBX) G4cout <<
	       Form("%4d: ID= %+5d trkId=%4d parentIdx=%4d (E, p)= %+7.2f/%+7.3f/%+7.3f/%+7.3f ekin= %11.8fMeV v=(%+8.3f, %+8.3f, %+14.8f), t = %f",
		    i, pdgid, tIdx, parentIdx, etot,
		    px, py, pz, ekin,
		    vx, vy, vz,
		    trj->GetGlobalTime()
		    )
		    << G4endl;

    // -- store production vertex of e+ from Mu (i.e. the Mu decay vertex)
    if ((11 == pdgid) && (1313 == TMath::Abs(getParticleID(parentIdx, evt)))) {
      TGenVtx *pVtx = rio->getEvent()->addGenVtx();
      pVtx->fNumber = rio->getEvent()->nGenVtx()-1;
      pVtx->fID = -131300;
      pVtx->fvMom.push_back(trj->GetParentID());
      pVtx->fvDau.push_back(trj->GetTrackID());
      pVtx->fLocalTime  = trj->GetLocalTime();
      pVtx->fGlobalTime = trj->GetGlobalTime();
      pVtx->fV.SetXYZ(trj->GetPoint(0)->GetPosition().x(),
		      trj->GetPoint(0)->GetPosition().y(),
		      trj->GetPoint(0)->GetPosition().z());
    }
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

// ----------------------------------------------------------------------
G4int EventAction::getParticleID(G4int idx, const G4Event *evt) {
  G4int mid(-1);

  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  for (G4int i = 0; i < n_trajectories; ++i) {
    //    G4Trajectory* trj=(G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
    Trajectory* trj=(Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
    if (idx == trj->GetTrackID()) {
      mid = trj->GetPDGEncoding();
      break;
    }
  }
  return mid;
}
