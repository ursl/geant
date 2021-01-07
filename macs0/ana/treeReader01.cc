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
  initVariables();

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
  ((TH1D*)fpHistFile->Get("h4"))->Fill(fpEvt->nGenCands());

  ((TH1D*)fpHistFile->Get("evts"))->Fill(0);
  if (nhtrk > 0) ((TH1D*)fpHistFile->Get("evts"))->Fill(2);
  if (nhmcp > 0) ((TH1D*)fpHistFile->Get("evts"))->Fill(3);
  if (nhtrk > 0 && nhmcp > 0) ((TH1D*)fpHistFile->Get("evts"))->Fill(4);

  TGenCand *pGen(0), *pEle(0), *pPos(0);
  double mass(-0.1);
  for (int igen = 0; igen < fpEvt->nGenCands(); ++igen) {
    pGen = fpEvt->getGenCand(igen);
    cout << "igen = " << igen << " pGen->fP.Vect()  = " << pGen->fP.Vect().X() << ","  << pGen->fP.Vect().Y() << ","  << pGen->fP.Vect().Z()
	 << " ptr: " << pGen
	 << endl;
    if (11 == pGen->fID) {
      pEle = pGen;
      fElePt = pEle->fP.Pt();
      fEleTheta = pEle->fP.Theta();
      fElePhi = pEle->fP.Phi();
      fEleE = pEle->fP.E();
      cout << "igen = " << igen << " pEle->fP.Vect()  = " << pEle->fP.X() << ","  << pEle->fP.Y() << ","  << pEle->fP.Z()
	   << " " << fElePt << "/" << fEleTheta << "/" << fElePhi
	   << " ptr: " << pEle
	   << endl;
    }
    if (-11 == pGen->fID) {
      pPos = pGen;
      fPosPt = pPos->fP.Pt();
      fPosTheta = pPos->fP.Theta();
      fPosPhi = pPos->fP.Phi();
      fPosE = pPos->fP.E();
      cout << "igen = " << igen << " pPos->fP.Vect()  = " << pPos->fP.X() << ","  << pPos->fP.Y() << ","  << pPos->fP.Z()
	   << " " << fPosPt << "/" << fPosTheta << "/" << fPosPhi
	   << " ptr: " << pPos
	   << endl;
      //      cout << "igen = " << igen << " pGen->fp.Pt() = " << pGen->fP.Pt() << endl;
    }
  }

  if (pEle && pPos) {
    TLorentzVector sum = pEle->fP + pPos->fP;
    mass = sum.M();
    fElePosMass = mass;
    fElePosOa = pPos->fP.Vect().Angle(pEle->fP.Vect());
    Double_t ptot2 = pPos->fP.Vect().Mag2()*pEle->fP.Vect().Mag2();
    Double_t arg1 = pPos->fP.Vect().Dot(pEle->fP.Vect());
    //    Double_t arg = arg1/TMath::Sqrt(ptot2);
    cout << "arg1 = " << arg1 << " ptot2 = " << ptot2 << endl
      //	 << " arg = " << arg
	 << " pEle->fP.Vect(1) = " << pEle->fP.X() << ","  << pEle->fP.Y() << ","  << pEle->fP.Z() << endl
	 << " pEle->fP.Vect(2) = " << pEle->fP.Vect().X() << ","  << pEle->fP.Vect().Y() << ","  << pEle->fP.Vect().Z() << endl
	 << " pPos->fP.Vect()  = " << pPos->fP.X() << ","  << pPos->fP.Y() << ","  << pPos->fP.Z() << endl
	 << " pPos->fP.Vect(2) = " << pPos->fP.Vect().X() << ","  << pPos->fP.Vect().Y() << ","  << pPos->fP.Vect().Z() << endl;
    if (pEle->fP.E() > 10 && pPos->fP.E() > 10) {
      ((TH1D*)fpHistFile->Get("mass"))->Fill(mass);
    }

    // ////////////////////////////////////////////////////////////////////////////////
    // /// Return the angle w.r.t. another 3-vector.

    // Double_t TVector3::Angle(const TVector3 & q) const
    // {
    //    Double_t ptot2 = Mag2()*q.Mag2();
    //    if(ptot2 <= 0) {
    //       return 0.0;
    //    } else {
    //       Double_t arg = Dot(q)/TMath::Sqrt(ptot2);
    //       if(arg >  1.0) arg =  1.0;
    //       if(arg < -1.0) arg = -1.0;
    //       return TMath::ACos(arg);
    //    }
    // }

    fpHistFile->cd();
    fillHist();
    fTree->Fill();
  }
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
  new TH1D("h4", "nGenCands", 40, 0., 40.);

  new TH1D("mass", "mass of e+e-", 40, -1., 99.);

  // -- Reduced Tree
  fTree = new TTree("events", "events");
  fTree->Branch("run",      &fRun,       "run/I");
  fTree->Branch("evt",      &fEvt,       "evt/I");

  fTree->Branch("elePt",    &fElePt,    "elePt/D");
  fTree->Branch("eleTheta", &fEleTheta, "eleTheta/D");
  fTree->Branch("elePhi",   &fElePhi,   "elePhi/D");
  fTree->Branch("eleE",     &fEleE,     "eleE/D");

  fTree->Branch("posPt",    &fPosPt,    "posPt/D");
  fTree->Branch("posTheta", &fPosTheta, "posTheta/D");
  fTree->Branch("posPhi",   &fPosPhi,   "posPhi/D");
  fTree->Branch("posE",     &fPosE,     "posE/D");

  fTree->Branch("eleposM",     &fElePosMass,"eleposM/D");
  fTree->Branch("eleposOa",    &fElePosOa,  "eleposOa/D");

}


// ----------------------------------------------------------------------
void treeReader01::initVariables() {
  cout << "treeReader01: initVariables: for run = " << fRun << "/evt = " << fEvt << endl;

  fElePt = -1.;
  fEleTheta = -1.;
  fElePhi = -1.;
  fEleE = -1.;
  fPosPt = -1.;
  fPosTheta = -1.;
  fPosPhi = -1.;
  fPosE = -1.;

  fElePosMass = fElePosOa = -1.;

}
