#include <iostream>

#include "rEvent.hh"


#define NGENCAND 100000

ClassImp(rEvent)

using namespace std;

// ----------------------------------------------------------------------
rEvent::rEvent() {
  fID = 0;

  fGenCands        = new TClonesArray("TGenCand", NGENCAND);
  fnGenCands       = 0;

}

// ----------------------------------------------------------------------
rEvent::~rEvent() {

}

// ----------------------------------------------------------------------
void rEvent::Clear(const char * /*opt*/) {
  clearGenBlock();
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

}


// ----------------------------------------------------------------------
TGenCand* rEvent::getGenCand(Int_t n) {
  return (TGenCand*)fGenCands->UncheckedAt(n);
}

// ----------------------------------------------------------------------
TGenCand* rEvent::addGenCand() {
  TClonesArray& d = *fGenCands;
  new(d[d.GetLast()+1]) TGenCand();
  ++fnGenCands;
  return (TGenCand*)d[d.GetLast()];
}

// ----------------------------------------------------------------------
void rEvent::dumpGenBlock() {
  TGenCand *pGenCand;
  for (int i = 0; i < fnGenCands; i++) {
    pGenCand = getGenCand(i);
    pGenCand->dump();
  }
}
