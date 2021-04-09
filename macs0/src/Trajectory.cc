#include "Trajectory.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4PrimaryParticle.hh"
#include "TrackInformation.hh"

#include "../common/musrMuonium.hh"



G4ThreadLocal G4Allocator<Trajectory> * myTrajectoryAllocator = 0;

// ----------------------------------------------------------------------
Trajectory::Trajectory(const G4Track* aTrack) : G4VTrajectory(), fPositionRecord(0), fParticleDefinition(0) {
  if (0) G4cout << "Trajectory::Trajectory!" << G4endl;
  fParticleDefinition = aTrack->GetDefinition();
  fParticleName = fParticleDefinition->GetParticleName();
  fPDGCharge = fParticleDefinition->GetPDGCharge();
  fPDGMass = fParticleDefinition->GetPDGMass();
  fPDGEncoding = fParticleDefinition->GetPDGEncoding();
  if ("unknown" == fParticleName) {
    G4PrimaryParticle*pp = aTrack->GetDynamicParticle()->GetPrimaryParticle();
    if (pp) {
      if(pp->GetCharge()<DBL_MAX) fPDGCharge = pp->GetCharge();
      fPDGEncoding = pp->GetPDGcode();
      if(pp->GetG4code()!=0) {
	fParticleName += " : ";
	fParticleName += pp->GetG4code()->GetParticleName();
      }
    }
  }
  fTrackID = aTrack->GetTrackID();
  TrackInformation* trackInfo = (TrackInformation*)(aTrack->GetUserInformation());
  fTrackStatus = trackInfo->GetTrackingStatus();
  if (2 == fTrackStatus) {
    fParentID = aTrack->GetParentID();
  } else {
    fParentID = -1;
  }
  fPositionRecord = new TrajectoryPointContainer();
  fPositionRecord->push_back(new G4TrajectoryPoint(aTrack->GetPosition()));
  fMomentum = aTrack->GetMomentum();
  fVertexPosition = aTrack->GetPosition();
  fGlobalTime = aTrack->GetGlobalTime();
  fLocalTime = aTrack->GetLocalTime();
}


// ----------------------------------------------------------------------
Trajectory::~Trajectory() {
  size_t i;
  for(i=0;i<fPositionRecord->size();i++){
    delete  (*fPositionRecord)[i];
  }
  fPositionRecord->clear();

  delete fPositionRecord;
}


// ----------------------------------------------------------------------
void Trajectory::ShowTrajectory(std::ostream& os) const
{
   os << G4endl << "TrackID =" << fTrackID
        << " : ParentID=" << fParentID << " : TrackStatus=" << fTrackStatus << G4endl;
   os << "Particle name : " << fParticleName << "  PDG code : " << fPDGEncoding
        << "  Charge : " << fPDGCharge << G4endl;
   os << "Original momentum : " <<
        G4BestUnit(fMomentum,"Energy") << G4endl;
   os << "Vertex : " << G4BestUnit(fVertexPosition,"Length")
        << "  Global time : " << G4BestUnit(fGlobalTime,"Time") << G4endl;
   os << "  Current trajectory has " << fPositionRecord->size()
        << " points." << G4endl;

   for( size_t i=0 ; i < fPositionRecord->size() ; i++){
       G4TrajectoryPoint* aTrajectoryPoint =
         (G4TrajectoryPoint*)((*fPositionRecord)[i]);
       os << "Point[" << i << "]"
            << " Position= " << aTrajectoryPoint->GetPosition() << G4endl;
   }
}

// ----------------------------------------------------------------------
void Trajectory::DrawTrajectory() const {

   G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
   G4ThreeVector pos;

   G4Polyline pPolyline;
   for (size_t i = 0; i < fPositionRecord->size() ; i++) {
     G4TrajectoryPoint* aTrajectoryPoint =
       (G4TrajectoryPoint*)((*fPositionRecord)[i]);
     pos = aTrajectoryPoint->GetPosition();
     pPolyline.push_back( pos );
   }

   G4Colour colour(0.2,0.2,0.2);
   if (G4Gamma::GammaDefinition() == fParticleDefinition)
     colour = G4Colour(0.99, 0.85, 0.0);  // gold
   else if (G4Electron::ElectronDefinition() == fParticleDefinition) {
     colour = G4Colour(0., 0., 0.8); // blue
   } else if (G4Positron::PositronDefinition() == fParticleDefinition) {
     colour = G4Colour(0.8, 0., 0.); // red
   } else if (G4MuonMinus::MuonMinusDefinition() == fParticleDefinition) {
     colour = G4Colour(0.99, 0.75, 0.85);   // pink
   } else if (G4MuonPlus::MuonPlusDefinition() == fParticleDefinition) {
     colour = G4Colour(0.99, 0.75, 0.85);  // pink
   } else if (musrMuonium::MuoniumDefinition() == fParticleDefinition) {
     colour = G4Colour(0., 1., 1.); // cyan
   } else if ("neutron" == fParticleDefinition->GetParticleType()) {
     colour = G4Colour(0.5, 0.5, 0.5);
   } else if ("neutrino" == fParticleDefinition->GetParticleType()) {
     colour = G4Colour(0.6, 0., 0.6);
   }

   G4VisAttributes attribs(colour);
   pPolyline.SetVisAttributes(attribs);
   if(pVVisManager) pVVisManager->Draw(pPolyline);
}

// ----------------------------------------------------------------------
const std::map<G4String,G4AttDef>* Trajectory::GetAttDefs() const {
  G4bool isNew;
  std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("Trajectory",isNew);
  if (isNew) {

    G4String id("ID");
    (*store)[id] = G4AttDef(id,"Track ID","Bookkeeping","","G4int");

    G4String pid("PID");
    (*store)[pid] = G4AttDef(pid,"Parent ID","Bookkeeping","","G4int");

    G4String status("Status");
    (*store)[status] = G4AttDef(status,"Track Status","Bookkeeping","","G4int");

    G4String pn("PN");
    (*store)[pn] = G4AttDef(pn,"Particle Name","Bookkeeping","","G4String");

    G4String ch("Ch");
    (*store)[ch] = G4AttDef(ch,"Charge","Physics","e+","G4double");

    G4String pdg("PDG");
    (*store)[pdg] = G4AttDef(pdg,"PDG Encoding","Bookkeeping","","G4int");

    G4String imom("IMom");
    (*store)[imom] = G4AttDef(imom, "Momentum of track at start of trajectory",
                              "Physics","G4BestUnit","G4ThreeVector");

    G4String imag("IMag");
    (*store)[imag] =
      G4AttDef(imag, "Magnitude of momentum of track at start of trajectory",
               "Physics","G4BestUnit","G4double");

    G4String vtxPos("VtxPos");
    (*store)[vtxPos] = G4AttDef(vtxPos, "Vertex position",
                              "Physics","G4BestUnit","G4ThreeVector");

    G4String ntp("NTP");
    (*store)[ntp] = G4AttDef(ntp,"No. of points","Bookkeeping","","G4int");

  }
  return store;
}

// ----------------------------------------------------------------------
std::vector<G4AttValue>* Trajectory::CreateAttValues() const {
  std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

  values->push_back
    (G4AttValue("ID",G4UIcommand::ConvertToString(fTrackID),""));

  values->push_back
    (G4AttValue("PID",G4UIcommand::ConvertToString(fParentID),""));

  values->push_back
    (G4AttValue("Status",G4UIcommand::ConvertToString(fTrackStatus),""));

  values->push_back(G4AttValue("PN",fParticleName,""));

  values->push_back
    (G4AttValue("Ch",G4UIcommand::ConvertToString(fPDGCharge),""));

  values->push_back
    (G4AttValue("PDG",G4UIcommand::ConvertToString(fPDGEncoding),""));

  values->push_back
    (G4AttValue("IMom",G4BestUnit(fMomentum,"Energy"),""));

  values->push_back
    (G4AttValue("IMag",G4BestUnit(fMomentum.mag(),"Energy"),""));

  values->push_back
    (G4AttValue("VtxPos",G4BestUnit(fVertexPosition,"Length"),""));

  values->push_back
    (G4AttValue("NTP",G4UIcommand::ConvertToString(GetPointEntries()),""));

  return values;
}

// ----------------------------------------------------------------------
void Trajectory::AppendStep(const G4Step* aStep) {
   fPositionRecord->push_back(
       new G4TrajectoryPoint(aStep->GetPostStepPoint()->GetPosition() ));
}

// ----------------------------------------------------------------------
void Trajectory::MergeTrajectory(G4VTrajectory* secondTrajectory) {
  if(!secondTrajectory) return;

  Trajectory* seco = (Trajectory*)secondTrajectory;
  G4int ent = seco->GetPointEntries();
  //
  // initial point of the second trajectory should not be merged
  for(int i=1;i<ent;i++)
  {
    fPositionRecord->push_back((*(seco->fPositionRecord))[i]);
  }
  delete (*seco->fPositionRecord)[0];
  seco->fPositionRecord->clear();

}
