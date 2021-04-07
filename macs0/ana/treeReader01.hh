#ifndef TREEREADER01_H
#define TREEREADER01_H

#include <iostream>
#include <string>

#include <TROOT.h>
#include <TString.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TTimeStamp.h>

#include "include/rEvent.hh"
#include "include/THit.hh"
#include "include/TGenCand.hh"
#include "include/TGenVtx.hh"

#define DR      57.29577951

class treeReader01 {
public:
  treeReader01(TChain *tree, TString evtClassName);
  virtual      ~treeReader01();
  virtual void       init(TString evtClassName);

  virtual void       openHistFile(TString filename);
  virtual void       closeHistFile();
  virtual void       bookHist();
  virtual void       readCuts(TString filename, int dump = 1);

  virtual void       startAnalysis();
  virtual void       endAnalysis();
  virtual int        loop(int nevents = 1, int start = -1);
  virtual TFile*     getFile() {return fpChain->GetCurrentFile();}
  virtual void       eventProcessing();
  virtual void       initVariables();
  virtual void       fillHist();
  virtual void       setVerbosity(int f) {std::cout << Form("setVerbosity(%d)", f) << std::endl;  fVerbose = f;}

  // -- study
  void fillMuFinal();
  void fillDaughters(TGenCand *pMu, int &idxEMuon, int &idxEAtom);

  // -- study in the context of mu3e for invariant mass of e+ Bhabha scattering on e-
  void doEnEpAnalysis();


  int fVerbose;

protected:

  TChain      *fpChain;        // pointer to the analyzed TTree or TChain
  TFile       *fpHistFile;     // for output histograms and reduced trees
  TString      fChainFileName; // the name of the chain file
  TString      fCutFile;       // contains file with the cut definitions

  rEvent *fpEvt;

  // -- Pre-filled variables
  int          fNentries;      // number of events in chain; filled in treeReader01::treeReader01()
  int          fChainEvent;    // current sequential event number in chain; filled in treeReader01::loop()
  int          fEvt;           // current event number; filled in treeReader01::loop()
  int          fRun;           // current run number; filled in treeReader01::loop()


  // -- vector with non-propagating Mu
  std::vector<TGenCand*> fMuFinal;

  // -- Histogram pointers
  TTree       *fTree;

  // -- Cut values
  double PTLO, PTHI;
  int TYPE;

  // -- variables
  double fElePt, fEleTheta, fElePhi, fEleE;
  double fPosPt, fPosTheta, fPosPhi, fPosE;

  double fElePosMass, fElePosOa;
};


#endif
