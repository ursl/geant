#ifndef INCLUDE_ROOTIO_HH
#define INCLUDE_ROOTIO_HH 1

#include <string>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

#include "TrackerHit.hh"
#include "MCPHit.hh"
#include "rEvent.hh"

class RootIO  {
public:
  virtual ~RootIO();

  static RootIO* GetInstance(std::string filename = "runG4.root");
  void WriteTrackerHits(std::vector<TrackerHit*>*);
  void WriteMCPHits(std::vector<MCPHit*>*);
  void Close();

  TTree* getTree() {return fTree;}
  rEvent* getEvent() {return fEvent;}
  void fillTree();
  void clear();

  void setVerbose(int v) {fVerbose = v;}

  static const int NGENMAX = 1000;


protected:
  RootIO(std::string filename = "g4run.root");

private:

  TFile*  fFile;
  TTree*  fTree;
  rEvent* fEvent;
  int fNevents;
  int fVerbose;

};
#endif // INCLUDE_ROOTIO_HH
