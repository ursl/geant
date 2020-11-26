#include "treeReader01.hh"

#include "TRandom.h"

using namespace std;

#include "treeReader01.icc"


// ----------------------------------------------------------------------
// Run with: ./runTreeReader01 -c chains/bg-test -D root
//           ./runTreeReader01 -f test.root
// ----------------------------------------------------------------------


// ----------------------------------------------------------------------
void treeReader01::startAnalysis() {
  cout << "treeReader01: startAnalysis: ..." << endl;
}

// ----------------------------------------------------------------------
void treeReader01::endAnalysis() {
  cout << "treeReader01: endAnalysis: ..." << endl;
}


// ----------------------------------------------------------------------
void treeReader01::eventProcessing() {

  cout << "----------------------------------------------------------------------" << endl;
  cout << "Found " << fpEvt->nGenCands() << " gen cands in event" << endl;
  cout << "Found " << fpEvt->nHits() << " hits in event" << endl;
  ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nHits());

  fpHistFile->cd();
  fillHist();

}


// ----------------------------------------------------------------------
void treeReader01::fillHist() {


}

// ----------------------------------------------------------------------
void treeReader01::bookHist() {
  cout << "==> treeReader01: bookHist> " << endl;

  new TH1D("h1", "nHits", 40, 0., 40.);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",    &fRun ,"run/I");

}


// ----------------------------------------------------------------------
void treeReader01::initVariables() {
  cout << "treeReader01: initVariables: ..." << endl;

  fRun = -1;
  fEvt = -1;
}
