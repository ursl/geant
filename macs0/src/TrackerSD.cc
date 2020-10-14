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
  collectionName.insert(name);
}

// ----------------------------------------------------------------------
TrackerSD::~TrackerSD() {
  //  RootIO::GetInstance()->Close();
}

// ----------------------------------------------------------------------
void TrackerSD::Initialize(G4HCofThisEvent* hce) {
  fHitsCollection = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  G4int hcID  = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );
}

// ----------------------------------------------------------------------
G4bool TrackerSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep==0.) return false;

  TrackerHit* newHit = new TrackerHit();
  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchable()
                                               ->GetReplicaNumber());
  newHit->SetEdep     (edep);
  newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());
  fHitsCollection->insert( newHit );

  return true;
}

// ----------------------------------------------------------------------
void TrackerSD::EndOfEvent(G4HCofThisEvent*) {
  // storing the hits in ROOT file
  G4int NbHits = fHitsCollection->entries();
  std::vector<TrackerHit*> hitsVector;

  {
    G4cout << "\n-------->Storing hits in the ROOT file: in this event there are " << NbHits
           << " hits in the tracker chambers: " << G4endl;
    for (G4int i=0;i<NbHits;i++) (*fHitsCollection)[i]->Print();
  }


  for (G4int i=0;i<NbHits;i++)
    hitsVector.push_back((*fHitsCollection)[i]);

  RootIO::GetInstance()->Write(&hitsVector);
}
