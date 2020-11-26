#include "THit.hh"
#include <iostream>

ClassImp(THit)

using namespace std;

THit::THit() { }

THit::THit(Int_t Option) { }

THit::THit(const  THit &other) {
  fNumber  = other.fNumber;
  fDetId   = other.fDetId;
  fChamber = other.fChamber;
  fTrack   = other.fTrack;
  fGenCand = other.fGenCand;
  fEdep    = other.fEdep;
  fPos     = other.fPos;
  fGblTime = other.fGblTime;
}


void THit::dump(int printV) {
  char line[200];
  if (1 == printV) {
    sprintf(line, "%4d %4d C%3d T%3d G%3d E=%8.3f vPos=(%+8.6f,%+8.6f,%+8.6f) t=%8.6f",
	    fNumber, fDetId, fChamber, fTrack, fGenCand,
	    fEdep, fPos.X(), fPos.Y(), fPos.Z(), fGblTime
	    );
  }
  cout << line << endl;
}

void THit::dump(ofstream &OUT, int printV) {
  char line[200];
  if (1 == printV) {
    sprintf(line, "%4d %4d C%3d T%3d G%3d E=%8.3f vPos=(%+8.6f,%+8.6f,%+8.6f) t=%8.6f",
	    fNumber, fDetId, fChamber, fTrack, fGenCand,
	    fEdep, fPos.X(), fPos.Y(), fPos.Z(), fGblTime
	    );
  }
  OUT << line << endl;

}
