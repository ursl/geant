#ifndef rEvent_h
#define rEvent_h 1

#include "TObject.h"
#include "TClonesArray.h"

#include "TGenCand.hh"
#include "TGenVtx.hh"
#include "THit.hh"


class rEvent: public TObject {
public:
  rEvent();
  ~rEvent();
  virtual void  Clear(Option_t *option ="");

  // ----------------------------------------------------------------------
  int                 nGenCands() {return fnGenCands;}
  TGenCand*           getGenCand(int n);
  TGenCand*           addGenCand();
  void                dumpGenBlock();
  void                clearGenBlock();

  // ----------------------------------------------------------------------
  int                 nGenVtx() {return fnGenVtx;}
  TGenVtx*            getGenVtx(int n);
  TGenVtx*            addGenVtx();
  void                dumpVtxBlock();
  void                clearVtxBlock();

  // ----------------------------------------------------------------------
  int                 nHits() {return fnHits;}
  THit*               getHit(int n);
  THit*               addHit();
  void                dumpHits();
  void                clearHits();

  // ----------------------------------------------------------------------
  // -- Basic event and detector information
  int                 fRunNumber, fEventNumber;
  // -- MC event/generation information
  int               fProcessID;

  // -- Reserve variables
  int               fIntRes[10];
  double            fDoubleRes[10];

  // ----------------------------------------------------------------------
  // -- the TClonesArray's
  int               fnGenCands;
  TClonesArray      *fGenCands;   //->

  int               fnGenVtx;
  TClonesArray      *fGenVtx;     //->

  int               fnHits;
  TClonesArray      *fHits;       //->

private:
  ClassDef(rEvent,1)


};

#endif
