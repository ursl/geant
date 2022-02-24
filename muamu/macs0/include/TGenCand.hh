#ifndef TGENCAND
#define TGENCAND


#include <fstream>

#include "TObject.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class TGenCand: public TObject {

public:

  TGenCand();
  TGenCand(int Option);
  TGenCand(const TGenCand &);
  ~TGenCand() { };
  void     clear() {fID = -123;}

  double  ekin() {return (fP.E() - fMass);}

  // ----------------------------------------------------------------------
  void dump(int mode = 2);
  void dump(std::ofstream &, int mode = 2);
  std::string dumpLine(int mode);

  // ----------------------------------------------------------------------
  int            fID, fNumber, fStatus;
  int            fMom1, fMom2;          // mothers
  int            fDau1, fDau2;          // daughters
  int            fTag;
  int 		 fQ;

  TLorentzVector fP;
  double         fMass;

  TVector3       fV;
  double         fLocalTime;
  double         fGlobalTime;

private:

  ClassDef(TGenCand,1)

};

#endif
