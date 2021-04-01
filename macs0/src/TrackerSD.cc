#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "RootIO.hh"

// ----------------------------------------------------------------------
TrackerSD::TrackerSD(const G4String& name, const G4String& hitsCollectionName) :
  G4VSensitiveDetector(name),
  fHitsCollection(NULL) {
  collectionName.insert(hitsCollectionName);
}

// ----------------------------------------------------------------------
TrackerSD::~TrackerSD() {
  //  RootIO::GetInstance()->Close();
}

// ----------------------------------------------------------------------
void TrackerSD::Initialize(G4HCofThisEvent* hce) {
  fHitsCollection = new TrackerHitsCollection(GetName(), collectionName[0]);
  G4int hcID  = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection(hcID, fHitsCollection);
}

// ----------------------------------------------------------------------
G4bool TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  G4double edep = aStep->GetTotalEnergyDeposit();

  if (edep==0.) return false;

  if (0) G4cout << "==========> TrackerSD::ProcessHits> new hit added " << G4endl;
  TrackerHit *newHit = new TrackerHit();
  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber());
  newHit->SetEdep     (edep);
  newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());
  newHit->SetEtrk     (aStep->GetTrack()->GetKineticEnergy());
  newHit->SetGblTime  (aStep->GetTrack()->GetGlobalTime());
  fHitsCollection->insert(newHit);

  return true;
}

// ----------------------------------------------------------------------
void TrackerSD::EndOfEvent(G4HCofThisEvent*) {
  // storing the hits in ROOT file
  G4int NbHits = fHitsCollection->entries();
  std::vector<TrackerHit*> hitsVector;

  G4cout << "-------->Storing hits in the ROOT file: there are " << NbHits
	 << " hits in the tracker chambers " << G4endl;
  if (0) {
    for (G4int i=0;i<NbHits;i++) {
      (*fHitsCollection)[i]->Print();
      // -- write into event/tree
      THit *hit = RootIO::GetInstance()->getEvent()->addHit();
      hit->fNumber  = RootIO::GetInstance()->getEvent()->nHits() - 1;
      hit->fDetId   = 0;// tracker = 0
      hit->fChamber = (*fHitsCollection)[i]->GetChamberNb();
      hit->fTrack   = (*fHitsCollection)[i]->GetTrackID();
      hit->fGenCand = 0;
      hit->fEdep    = (*fHitsCollection)[i]->GetEdep();
      hit->fGblTime = (*fHitsCollection)[i]->GetGblTime();
      G4ThreeVector x = (*fHitsCollection)[i]->GetPos();
      hit->fPos     = TVector3(x.x(), x.y(), x.z());
    }
  }

  if (0) {
    for (G4int i=0;i<NbHits;i++) hitsVector.push_back((*fHitsCollection)[i]);
    RootIO::GetInstance()->WriteTrackerHits(&hitsVector);
  }
}
