#ifndef rEvent_h
#define rEvent_h 1

#include "TObject.h"
#include "TClonesArray.h"

#include "TGenCand.hh"


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


  int fID;

  int               fnGenCands;
  TClonesArray      *fGenCands;         //->

private:
  ClassDef(rEvent,1)


};

#endif
