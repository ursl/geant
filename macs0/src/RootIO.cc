#include <sstream>

#include "RootIO.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

static RootIO* instance = 0;

using namespace std;

// ----------------------------------------------------------------------
RootIO::RootIO(string filename):fNevents(0), fVerbose(0) {
  gSystem->Load("libClassesDict");

  fFile = TFile::Open(filename.c_str(), "RECREATE");


  fTree = new TTree("gen", "gen");
  fTree->SetDirectory(fFile);

  fEvent = new rEvent;
  rEvent::Class()->SetCanSplit(1);
  fTree->Branch("rEvent", "rEvent", &fEvent, 64000, 99);
}

// ----------------------------------------------------------------------
RootIO::~RootIO() {
  Close();
}


// ----------------------------------------------------------------------
RootIO* RootIO::GetInstance(string filename) {
  if (instance == 0)   {
    instance = new RootIO(filename);
  }
  return instance;
}

// ----------------------------------------------------------------------
void RootIO::clear() {
}

// ----------------------------------------------------------------------
void RootIO::fillTree() {
  if (fVerbose > 0) G4cout << "Filling tree with ngen = " << fEvent->nGenCands() << G4endl;
  if (0) {
    fTree->Print();
  }
  fTree->Fill();
}

// ----------------------------------------------------------------------
void RootIO::WriteTrackerHits(std::vector<TrackerHit*>* hcont) {
  fNevents++;

  std::ostringstream os;
  os << fNevents;
  std::string stevt = "Event_" + os.str();
  const char* chevt = stevt.c_str();

  if (fVerbose > 0) G4cout << "writing " << stevt << G4endl;
  fFile->WriteObject(hcont, chevt);
}

// ----------------------------------------------------------------------
void RootIO::WriteMCPHits(std::vector<MCPHit*>* hcont) {
  fNevents++;

  std::ostringstream os;
  os << fNevents;
  std::string stevt = "Event_" + os.str();
  const char* chevt = stevt.c_str();

  if (fVerbose > 0) G4cout << "writing " << stevt << G4endl;
  fFile->WriteObject(hcont, chevt);
}

// ----------------------------------------------------------------------
void RootIO::Close() {
  G4cout << "RootIO::Close() " << G4endl;

  fTree->Write();
  fFile->Write();

  fFile->Close();
}
