#ifndef rEvent_h
#define rEvent_h 1

#include "TObject.h"

class rEvent: public TObject {
public:
  rEvent();
  ~rEvent();

  int fID;
};

#endif
