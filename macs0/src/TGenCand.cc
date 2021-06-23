#include "TGenCand.hh"
#include <iostream>

ClassImp(TGenCand)

using namespace std;

TGenCand::TGenCand() { }

TGenCand::TGenCand(Int_t Option) { }

TGenCand::TGenCand(const  TGenCand &other) {
  fID         = other.fID;
  fNumber     = other.fNumber;
  fStatus     = other.fStatus;
  fMom1       = other.fMom1;
  fMom2       = other.fMom2;
  fDau1       = other.fDau1;
  fDau2       = other.fDau2;
  fTag        = other.fTag;
  fQ          = other.fQ;
  fP          = other.fP;
  fMass       = other.fMass;
  fV          = other.fV;
  fGlobalTime = other.fGlobalTime;
  fLocalTime  = other.fLocalTime;
}


// ----------------------------------------------------------------------
string TGenCand::dumpLine(int mode) {
  char line[200];
  if (2 == mode) {
    // -- Ekin(px, py, pz)
    sprintf(line, "%4d %+6d S%2d mom(%4d,%4d) dau(%5d,%5d) Ekin=%10.7f(%+9.3f,%+9.3f,%+9.3f) v=(%+10.4f,%+10.4f,%+13.6f), t=%+16.11f",
	    fNumber, fID, fStatus, fMom1, fMom2, fDau1, fDau2,
	    fP.E() - fMass,
	    fP.X(), fP.Y(), fP.Z(),
	    fV.X(), fV.Y(), fV.Z(), fGlobalTime);
  } else if (1 == mode) {
    // -- pT(pT, eta, phi)
    sprintf(line, "%4d %+6d S%2d mom(%4d,%4d) dau(%5d,%5d) p/t=%8.3f(%+9.3f,%+9.3f,%+9.3f) v=(%+10.4f,%+10.4f,%+13.6f), t=%+16.11f",
	    fNumber, fID, fStatus, fMom1, fMom2, fDau1, fDau2,
	    fP.Rho(),
	    fP.Perp(), (fP.Perp()>0.01?fP.Eta():99.), fP.Phi(),
	    fV.X(), fV.Y(), fV.Z(), fGlobalTime);
  } else {
    // -- p(px, py, pz)
    sprintf(line, "%4d %+6d S%2d mom(%4d,%4d) dau(%5d,%5d) p/c=%8.3f(%+9.3f,%+9.3f,%+9.3f) v=(%+10.4f,%+10.4f,%+13.6f), t=%+16.11f",
	    fNumber, fID, fStatus, fMom1, fMom2, fDau1, fDau2,
	    fP.Rho(),
	    fP.X(), fP.Y(), fP.Z(),
	    fV.X(), fV.Y(), fV.Z(), fGlobalTime);

  }
  return string(line);
}

// ----------------------------------------------------------------------
void TGenCand::dump(int mode) {
  cout << dumpLine(mode) << endl;
}

void TGenCand::dump(ofstream &OUT, int mode) {
  OUT << dumpLine(mode) << endl;
}
