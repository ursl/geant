#ifndef INCLUDE_ROOTIO_HH
#define INCLUDE_ROOTIO_HH 1

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

  static RootIO* GetInstance();
  void WriteTrackerHits(std::vector<TrackerHit*>*);
  void WriteMCPHits(std::vector<MCPHit*>*);
  void Close();

  TTree* getTree() {return fTree;}
  rEvent* getEvent() {return fEvent;}
  void fillTree();
  void clear();

  static const int NGENMAX = 100;


protected:
  RootIO();

private:

  TFile*  fFile;
  TTree*  fTree;
  rEvent* fEvent;
  int fNevents;


};
#endif // INCLUDE_ROOTIO_HH
