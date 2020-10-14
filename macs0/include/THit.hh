#ifndef THIT
#define THIT


#include <fstream>

#include "TObject.h"
#include "TVector3.h"

class THit: public TObject {

public:

  THit();
  THit(int Option);
  THit(const THit &);
  ~THit() { };
  void     clear() {fID = -123;}

  // ----------------------------------------------------------------------
  void dump(int printV = 1);
  void dump(std::ofstream &, int printV = 1);

  // ----------------------------------------------------------------------
  int            fNumber, fID, fChamber, fTrack, fGenCand;
  double         fEdep;
  TVector3       fPos;
  double         fTime;

private:

  ClassDef(THit,1)

};

#endif
