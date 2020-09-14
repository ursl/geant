#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "RootIO.hh"

// ----------------------------------------------------------------------
TrackerSD::TrackerSD(G4String name) :G4VSensitiveDetector(name), fTrackerCollection(0), fHCID(0) {
  G4String HCname = name + "_HC";
  collectionName.insert(HCname);
  G4cout << collectionName.size() << "   CalorimeterSD name:  " << name << " collection Name: "
         << HCname << G4endl;
  fHCID = -1;
}

// ----------------------------------------------------------------------
TrackerSD::~TrackerSD() {
  RootIO::GetInstance()->Close();
}

// ----------------------------------------------------------------------
void TrackerSD::Initialize(G4HCofThisEvent* HCE) {
  fTrackerCollection = new TrackerHitsCollection
    (SensitiveDetectorName,collectionName[0]);
  if (fHCID < 0) {
    G4cout << "CalorimeterSD::Initialize:  " << SensitiveDetectorName << "   "
           << collectionName[0] << G4endl;
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);

  }
  HCE->AddHitsCollection(fHCID, fTrackerCollection);
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
  fTrackerCollection->insert( newHit );

  return true;
}

// ----------------------------------------------------------------------
void TrackerSD::EndOfEvent(G4HCofThisEvent*) {
  // storing the hits in ROOT file
  G4int NbHits = fTrackerCollection->entries();
  std::vector<TrackerHit*> hitsVector;

  {
    G4cout << "\n-------->Storing hits in the ROOT file: in this event there are " << NbHits
           << " hits in the tracker chambers: " << G4endl;
    for (G4int i=0;i<NbHits;i++) (*fTrackerCollection)[i]->Print();
  }


  for (G4int i=0;i<NbHits;i++)
    hitsVector.push_back((*fTrackerCollection)[i]);

  RootIO::GetInstance()->Write(&hitsVector);
}
