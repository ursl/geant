#include "TGenVtx.hh"
#include <iostream>

ClassImp(TGenVtx)

using namespace std;

TGenVtx::TGenVtx() { }

TGenVtx::TGenVtx(Int_t Option) { }

TGenVtx::TGenVtx(const  TGenVtx &other) {
  fID     = other.fID;
  fNumber = other.fNumber;
  fStatus = other.fStatus;
  fTag    = other.fTag;
  fQ      = other.fQ;
  fV      = other.fV;
  fLocalTime  = other.fLocalTime;
  fGlobalTime = other.fGlobalTime;
}

TGenVtx::~TGenVtx() {
  clear();
}

// ----------------------------------------------------------------------
void TGenVtx::clear() {
  int vinit[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
  fvMom = std::vector<int>(vinit, vinit + sizeof(vinit) / sizeof(int) );
  fvDau = std::vector<int>(vinit, vinit + sizeof(vinit) / sizeof(int) );

  fID = fNumber = fStatus = fTag = -1;
  fQ = -9999;
  fLocalTime = -9999.;
  fGlobalTime = -9999.;
  fV.SetXYZ(-1e30,-1e30,-1e30);
}


// ----------------------------------------------------------------------
void TGenVtx::dump(int printOption) {
  char line[200];
  if (1 == printOption) {
    sprintf(line, "%4d %+6d mom=%3d v=(%+10.4f,%+10.4f,%+13.6f) t = %10.8f/%10.8f",
	    fNumber, fID, fvMom.at(0), fV.X(), fV.Y(), fV.Z(), fGlobalTime, fLocalTime);
    cout << line << endl;
  }
}

void TGenVtx::dump(ofstream &OUT) {
  char line[200];
  sprintf(line, "%4d %+6d v=(%+10.4f,%+10.4f,%+13.6f) t = %10.8f/%10.8f",
	  fNumber, fID, fV.X(), fV.Y(), fV.Z(), fGlobalTime, fLocalTime);
  OUT << line << endl;
}
