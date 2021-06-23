#include <iostream>

#include "rEvent.hh"


#define NGENCAND 10000
#define NGENVTX  1000
#define NHITS 10000

ClassImp(rEvent)

using namespace std;

// ----------------------------------------------------------------------
rEvent::rEvent() {

  fGenCands  = new TClonesArray("TGenCand", NGENCAND);
  fGenVtx    = new TClonesArray("TGenVtx", NGENVTX);
  fHits      = new TClonesArray("THit", NHITS);

  Clear();

}

// ----------------------------------------------------------------------
rEvent::~rEvent() {

}

// ----------------------------------------------------------------------
void rEvent::Clear(const char * /*opt*/) {
  for (int i = 0; i < 10; ++i) {
    fIntRes[i]    = 0;
    fDoubleRes[i] = 0.;
  }

  fRunNumber = fEventNumber = fProcessID = 0;

  clearGenBlock();
  clearHits();
}

// ----------------------------------------------------------------------
void rEvent::clearGenBlock() {
  TGenCand *pGenCand;
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    pGenCand->clear();
  }
  fGenCands->Clear();
  fnGenCands = 0;

  TGenVtx *pGenVtx;
  for (int i = 0; i < fnGenVtx; i++) {
    pGenVtx = getGenVtx(i);
    pGenVtx->clear();
  }
  fGenVtx->Clear();
  fnGenVtx = 0;

}


// ----------------------------------------------------------------------
TGenCand* rEvent::getGenCand(Int_t n) {
  return (TGenCand*)fGenCands->UncheckedAt(n);
}


// ----------------------------------------------------------------------
TGenCand* rEvent::getGenCandWithNumber(int number) {
  TGenCand *pGenCand;
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    if (number == pGenCand->fNumber) {
      return (TGenCand*)fGenCands->UncheckedAt(i);
    }
  }
  return 0;
}


// ----------------------------------------------------------------------
TGenVtx* rEvent::getGenVtx(Int_t n) {
  return (TGenVtx*)fGenVtx->UncheckedAt(n);
}


// ----------------------------------------------------------------------
TGenCand* rEvent::addGenCand() {
  TClonesArray& d = *fGenCands;
  new(d[d.GetLast()+1]) TGenCand();
  ++fnGenCands;
  return (TGenCand*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
TGenVtx* rEvent::addGenVtx() {
  TClonesArray& d = *fGenVtx;
  new(d[d.GetLast()+1]) TGenVtx();
  ++fnGenVtx;
  return (TGenVtx*)d[d.GetLast()];
}


// ----------------------------------------------------------------------
void rEvent::dumpGenBlock() {
  TGenCand *pGenCand;
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    pGenCand->dump();
  }

  TGenVtx *pGenVtx;
  for (int i = 0; i < fnGenVtx; i++) {
    pGenVtx = getGenVtx(i);
    pGenVtx->dump();
  }
}


// ----------------------------------------------------------------------
void rEvent::clearHits() {
  THit *pHit;
  for (int i = 0; i < fnHits; i++) {
    pHit = getHit(i);
    pHit->clear();
  }
  fHits->Clear();
  fnHits = 0;
}

// ----------------------------------------------------------------------
THit* rEvent::getHit(Int_t n) {
  return (THit*)fHits->UncheckedAt(n);
}

// ----------------------------------------------------------------------
THit* rEvent::addHit() {
  TClonesArray& d = *fHits;
  new(d[d.GetLast()+1]) THit();
  ++fnHits;
  return (THit*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
void rEvent::dumpHits() {
  THit *pHit;
  for (int i = 0; i < fnHits; i++) {
    pHit = getHit(i);
    pHit->dump();
  }
}
