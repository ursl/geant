#include "treeReader01.hh"

#include "TRandom.h"

#define MMUON 105.658305

using namespace std;


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

  // -- generic rudimentary analysis
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


  // -- call specific analysis functions
  TGenCand *pGen(0);
  for (int igen = 0; igen < fpEvt->nGenCands(); ++igen) {
    pGen = fpEvt->getGenCand(igen);
    if ((0 == pGen->fMom1) && (-13 == pGen->fID)) {
      double gamma = pGen->fP.E()/MMUON;
      //      double ekin = (gamma-1)*MMUON;
      double ekin = 1.e3*(pGen->fP.E() - MMUON); // in keV
      cout << "gamma = " << gamma << " ekin = " << ekin << endl;
      ((TH1D*)fpHistFile->Get("h7"))->Fill(ekin);
    }
    pGen->dump(0);
  }
  fillMuFinal();
  fillHist();
  cout << "----------------------------------------------------------------------" << endl;
  cout << "vtx: " << endl;
  TGenVtx *pVtx;
  for (int ivtx = 0; ivtx < fpEvt->nGenVtx(); ++ivtx) {
    pVtx = fpEvt->getGenVtx(ivtx);
    pVtx->dump();
  }
}

// ----------------------------------------------------------------------
void treeReader01::fillMuFinal() {
  fMuFinal.clear();

  vector<TGenCand*> mus;
  TGenCand *pGen(0), *pMu(0), *pEle(0), *pPos(0), *pMom(0), *pDau(0);
  double mass(-0.1);
  for (int igen = 0; igen < fpEvt->nGenCands(); ++igen) {
    pGen = fpEvt->getGenCand(igen);
    if (-1313 == pGen->fID) {
      // -- undecayed Mu (propagating Mu have one daughter and therefore fDau1 not unfilled)
      if (9999 == pGen->fDau1  && -9999 == pGen->fDau2) {
	fMuFinal.push_back(pGen);
      }
      // -- Mu that are not undecayed and have more than 1 daughter
      if ((pGen->fDau2 - pGen->fDau1 > 0) && (pGen->fDau1 < 9999) && (pGen->fDau2 > -9999)) {
	fMuFinal.push_back(pGen);
      }
    }
  }

  if (1) {
    cout << "and the final Mu list" << endl;
    for (unsigned int i = 0; i < fMuFinal.size(); ++i) {
      fMuFinal[i]->dump(0);
    }
  }
}

// ----------------------------------------------------------------------
void treeReader01::fillDaughters(TGenCand *pMu, int &idxEMuon, int &idxEAtom) {
  if (9999 == pMu->fDau1) return;
  if (-9999 == pMu->fDau2) return;

  cout << "fillDaughters: " << idxEAtom << "  " << idxEMuon << endl;
  // -- the numerical index does not at all correspond to the trkID (=fNumber)
  TGenCand *pGen(0);
  vector<int> vDau;
  for (unsigned int j = pMu->fDau1; j <= pMu->fDau2; ++j) {
    for (int i = 0; i < fpEvt->nGenCands(); ++i) {
      pGen = fpEvt->getGenCand(i);
      if (j == pGen->fNumber) {
	vDau.push_back(i);
	break;
      }
    }
  }

  for (unsigned int i = 0; i < vDau.size(); ++i) {
    cout << "look at vDau[i]->fNumber = " << fpEvt->getGenCand(vDau[i])->fNumber
	 << " vDau[i]->fID = " << fpEvt->getGenCand(vDau[i])->fID
	 << endl;
    if (-11 == fpEvt->getGenCand(vDau[i])->fID) {
      idxEMuon = vDau[i];
      cout << "  setting pEMuon = " << fpEvt->getGenCand(vDau[i])
	   << "  idxEMuon = " << idxEMuon
	   << endl;
    }
    if (+11 == fpEvt->getGenCand(vDau[i])->fID) {
      idxEAtom = vDau[i];
      cout << "  setting pEAtom = " << fpEvt->getGenCand(vDau[i])
	   << "  idxEAtom = " << idxEAtom
	   << endl;
    }
  }
  cout  << "pMu->fNumber = " << pMu->fNumber
	<< "  idxEMuon = " << idxEMuon
	<< "  idxEAtom = " << idxEAtom
	<< endl;
}

// ----------------------------------------------------------------------
void treeReader01::fillHist() {
  int idxEAtom(-1), idxEMuon(-1);
  TGenCand *pMu(0), *pEMuon(0), *pEAtom(0);
  for (unsigned int i = 0; i < fMuFinal.size(); ++i) {
    pMu = fMuFinal[i];
    fillDaughters(pMu, idxEMuon, idxEAtom);
    if (idxEMuon > -1) pEMuon = fpEvt->getGenCand(idxEMuon);
    if (idxEAtom > -1) pEAtom = fpEvt->getGenCand(idxEAtom);
    cout << "Filling MuFinal[" << i << "] at " << fMuFinal[i]->fNumber
	 << " with idxEMuon = " << idxEMuon
	 << " with idxEAtom = " << idxEAtom
	 << " pEMuon = " << pEMuon << " pEAtom = " << pEAtom
	 << " pMu->fMass = " << pMu->fMass
	 << endl;
    if (0 == pEMuon) continue;
    if (0 == pEAtom) continue;
    if (pEMuon) {
      double zdecay = pEMuon->fV.Z();
      double fl = (pMu->fV - pEMuon->fV).Mag();
      double t  = MMUON*fl/pMu->fP.Rho()/(1000.*TMath::C()); // mm
      ((TH1D*)fpHistFile->Get("h5"))->Fill(zdecay);
      ((TH1D*)fpHistFile->Get("h8"))->Fill(fl);
      ((TH1D*)fpHistFile->Get("h9"))->Fill(pMu->fP.Rho());
      ((TH2D*)fpHistFile->Get("m10"))->Fill(t, zdecay);
      if (zdecay > -400 && zdecay < 2095.) {
	((TH1D*)fpHistFile->Get("h10"))->Fill(t);
      }
    }
    ((TH1D*)fpHistFile->Get("h6"))->Fill(pMu->fV.Z());
  }
}

// ----------------------------------------------------------------------
void treeReader01::bookHist() {
  cout << "==> treeReader01: bookHist> " << endl;

  new TH1D("evts", "events", 40, 0., 40.);
  new TH1D("h1", "nHits", 40, 0., 40.);
  new TH1D("h2", "nHits Trk", 40, 0., 40.);
  new TH1D("h3", "nHits MCP", 40, 0., 40.);
  new TH1D("h4", "nGenCands", 40, 0., 40.);
  new TH1D("h5", "z(Mu decay) [mm]", 500, -2200., 2200.);
  new TH1D("h6", "z(Mu produced) [mm]", 500, -2200., 2200.);

  new TH1D("h7", "mu+(beam) Ekin [keV]", 100, 0., 100.);
  new TH1D("h8", "Mu decay length [mm]", 100, 0., 5000.);
  new TH1D("h9", "Mu momentum [MeV]", 100, 0., 10.);
  new TH1D("h10", "proper decay time", 100, 0., 1.e-5);
  new TH2D("m10", "proper decay time vs. z", 100, 0., 1.e-5, 100, -2200., 2200.);

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


// ----------------------------------------------------------------------
void treeReader01::doEnEpAnalysis() {
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




// ======================================================================
// -- Below is the icc material
// ======================================================================

// ----------------------------------------------------------------------
treeReader01::treeReader01(TChain *chain, TString evtClassName) {
  cout << "==> treeReader01: constructor..." << endl;
  if (chain == 0) {
    cout << "You need to pass a chain!" << endl;
  }
  fpChain = chain;
  fNentries = chain->GetEntries();
  init(evtClassName);
}

// ----------------------------------------------------------------------
void treeReader01::init(TString evtClassName) {
  fpEvt = new rEvent();
  cout << "==> treeReader01: init ..." << fpEvt << endl;
  fpChain->SetBranchAddress(evtClassName, &fpEvt);
  initVariables();
}

// ----------------------------------------------------------------------
treeReader01::~treeReader01() {
  cout << "==> treeReader01: destructor ..." << endl;
  if (!fpChain) return;
  delete fpChain->GetCurrentFile();
}

// ----------------------------------------------------------------------
void treeReader01::openHistFile(TString filename) {
  fpHistFile = new TFile(filename.Data(), "RECREATE");
  fpHistFile->cd();
  cout << "==> treeReader01: Opened " << fpHistFile->GetName() << endl;
}

// ----------------------------------------------------------------------
void treeReader01::closeHistFile() {
  cout << "==> treeReader01: Writing " << fpHistFile->GetName() << endl;
  fpHistFile->cd();
  fpHistFile->Write();
  fpHistFile->Close();
  delete fpHistFile;

}

// --------------------------------------------------------------------------------------------------
void treeReader01::readCuts(TString filename, int dump) {
  char  buffer[200];
  fCutFile = filename;
  if (dump) cout << "==> treeReader01: Reading " << fCutFile.Data() << " for cut settings" << endl;
  sprintf(buffer, "%s", fCutFile.Data());
  ifstream is(buffer);
  char CutName[100];
  float CutValue;
  int ok(0);

  TString fn(fCutFile.Data());

  if (dump) {
    cout << "====================================" << endl;
    cout << "==> treeReader01: Cut file  " << fCutFile.Data() << endl;
    cout << "------------------------------------" << endl;
  }

  TH1D *hcuts = new TH1D("hcuts", "", 1000, 0., 1000.);
  hcuts->GetXaxis()->SetBinLabel(1, fn.Data());
  int ibin;

  while (is.getline(buffer, 200, '\n')) {
    ok = 0;
    if (buffer[0] == '#') {continue;}
    if (buffer[0] == '/') {continue;}
    sscanf(buffer, "%s %f", CutName, &CutValue);

    if (!strcmp(CutName, "TYPE")) {
      TYPE = int(CutValue); ok = 1;
      if (dump) cout << "TYPE:           " << TYPE << endl;
    }

    if (!strcmp(CutName, "PTLO")) {
      PTLO = CutValue; ok = 1;
      if (dump) cout << "PTLO:           " << PTLO << " GeV" << endl;
      ibin = 11;
      hcuts->SetBinContent(ibin, PTLO);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{min}(l) [GeV]");
    }

    if (!strcmp(CutName, "PTHI")) {
      PTHI = CutValue; ok = 1;
      if (dump) cout << "PTHI:           " << PTHI << " GeV" << endl;
      ibin = 12;
      hcuts->SetBinContent(ibin, PTHI);
      hcuts->GetXaxis()->SetBinLabel(ibin, "p_{T}^{max}(l) [GeV]");
    }


    if (!ok) cout << "==> treeReader01: ERROR: Don't know about variable " << CutName << endl;
  }

  if (dump)  cout << "------------------------------------" << endl;
}


// ----------------------------------------------------------------------
int treeReader01::loop(int nevents, int start) {
  int nb = 0, maxEvents(0);

  cout << "==> treeReader01: Chain has a total of " << fNentries << " events" << endl;

  // -- Setup for restricted running (not yet foolproof, i.e. bugfree)
  if (nevents < 0) {
    maxEvents = fNentries;
  } else {
    cout << "==> treeReader01: Running over " << nevents << " events" << endl;
    maxEvents = nevents;
  }
  if (start < 0) {
    start = 0;
  } else {
    cout << "==> treeReader01: Starting at event " << start << endl;
    if (maxEvents >  fNentries) {
      cout << "==> treeReader01: Requested to run until event " << maxEvents << ", but will run only to end of chain at ";
      maxEvents = fNentries;
      cout << maxEvents << endl;
    } else {
      cout << "==> treeReader01: Requested to run until event " << maxEvents << endl;
    }
  }

  // -- The main loop
  int step(50000);
  if (maxEvents < 1000000) step = 10000;
  if (maxEvents < 100000)  step = 1000;
  if (maxEvents < 10000)   step = 500;
  if (maxEvents < 1000)    step = 100;

  int treeNumber(0), oldTreeNumber(-1);
  fpChain->GetFile(); // without this, treeNumber initially will be -1.
  for (int jEvent = start; jEvent < maxEvents; ++jEvent) {
    treeNumber = fpChain->GetTreeNumber();
    if (treeNumber != oldTreeNumber) {
      cout << "    " << Form("      %8d", jEvent) << " " << fpChain->GetFile()->GetName() << endl;
      oldTreeNumber = treeNumber;
    }

    if (jEvent%step == 0) cout << Form(" .. Event %8d", jEvent);

    fChainEvent = jEvent;
    fpEvt->Clear();
    nb += fpChain->GetEvent(jEvent);

    fEvt = static_cast<long int>(fpEvt->fEventNumber);
    fRun = static_cast<long int>(fpEvt->fRunNumber);

    if (jEvent%step == 0) {
      TTimeStamp ts;
      cout  << " (run: " << Form("%8d", fRun)
	    << ", event: " << Form("%10d", fEvt)
	    << ", time now: " << ts.AsString("lc")
	    << ")" << endl;
    }

    eventProcessing();
  }
  return 0;

}
