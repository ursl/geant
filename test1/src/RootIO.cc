#include <sstream>

#include "RootIO.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

static RootIO* instance = 0;

// ----------------------------------------------------------------------
RootIO::RootIO() {
  fNevents = 0;
  TSystem ts;
  gSystem->Load("libClassesDict");

  fFile = TFile::Open("g4run3.root","RECREATE");


  fTree = new TTree("gen", "gen");
  fTree->SetDirectory(fFile);

  fEvent = new rEvent;

  fTree->Branch("event", &fEvent);
}

// ----------------------------------------------------------------------
RootIO::~RootIO() {
  Close();
}

// ----------------------------------------------------------------------
RootIO* RootIO::GetInstance() {
  if (instance == 0)   {
    instance = new RootIO();
  }
  return instance;
}

// ----------------------------------------------------------------------
void RootIO::clear() {
  fNGenParticles = -99;
}

// ----------------------------------------------------------------------
void RootIO::fillTree() {
  G4cout << "Filling tree with ngen = " << fNGenParticles << G4endl;
  fTree->Print();
  fTree->Fill();
}

// ----------------------------------------------------------------------
void RootIO::Write(std::vector<TrackerHit*>* hcont) {
  fNevents++;

  std::ostringstream os;
  os << fNevents;
  std::string stevt = "Event_" + os.str();
  const char* chevt = stevt.c_str();

  G4cout << "writing " << stevt << G4endl;
  fFile->WriteObject(hcont, chevt);
}

// ----------------------------------------------------------------------
void RootIO::Close() {
  G4cout << "RootIO::Close() " << G4endl;

  fTree->Write();
  fFile->Write();

  fFile->Close();
}
