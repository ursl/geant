#ifndef TGENVTX
#define TGENVTX


#include <fstream>

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include <vector>

class TGenVtx: public TObject {

public:

  TGenVtx();
  TGenVtx(int Option);
  TGenVtx(const TGenVtx &);
  ~TGenVtx();
  void     clear();

  // ----------------------------------------------------------------------
  void dump(int printPt = 1);
  void dump(std::ofstream &);

  // ----------------------------------------------------------------------
  int            fNumber, fID, fStatus, fTag, fQ;
  std::vector<int> fvMom;          // mother indices
  std::vector<int> fvDau;          // daughter indices

  TVector3       fV;
  double         fLocalTime, fGlobalTime;

private:

  ClassDef(TGenVtx,1)

};

#endif
