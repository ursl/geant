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

  if (0) {
    cout << "----------------------------------------------------------------------" << endl;
    cout << "Found " << fpEvt->nGenCands() << " gen cands in event" << endl;
    cout << "Found " << fpEvt->nHits() << " hits in event" << endl;
  }
  ((TH1D*)fpHistFile->Get("h1"))->Fill(fpEvt->nHits());
  THit *pHit(0);
  int nhtrk(0), nhmcp(0);
  for (int ihit = 0; ihit < fpEvt->nHits(); ++ihit) {
    pHit = fpEvt->getHit(ihit);
    if (0 == pHit->fDetId) ++nhtrk;
    if (1 == pHit->fDetId) ++nhmcp;
  }
  ((TH1D*)fpHistFile->Get("h2"))->Fill(nhtrk);
  ((TH1D*)fpHistFile->Get("h3"))->Fill(nhmcp);

  ((TH1D*)fpHistFile->Get("evts"))->Fill(0);
  if (nhtrk > 0) ((TH1D*)fpHistFile->Get("evts"))->Fill(2);
  if (nhmcp > 0) ((TH1D*)fpHistFile->Get("evts"))->Fill(3);
  if (nhtrk > 0 && nhmcp > 0) ((TH1D*)fpHistFile->Get("evts"))->Fill(4);

  fpHistFile->cd();
  fillHist();

}


// ----------------------------------------------------------------------
void treeReader01::fillHist() {


}

// ----------------------------------------------------------------------
void treeReader01::bookHist() {
  cout << "==> treeReader01: bookHist> " << endl;

  new TH1D("evts", "events", 40, 0., 40.);
  new TH1D("h1", "nHits", 40, 0., 40.);
  new TH1D("h2", "nHits Trk", 40, 0., 40.);
  new TH1D("h3", "nHits MCP", 40, 0., 40.);

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
