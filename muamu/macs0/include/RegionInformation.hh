#ifndef RegionInformation_H
#define RegionInformation_H 1

#include "globals.hh"
#include "G4VUserRegionInformation.hh"

class RegionInformation : public G4VUserRegionInformation
{
public:
  RegionInformation();
  virtual ~RegionInformation();
  virtual void Print() const;

  inline void SetWorld(G4bool v=true) {fIsWorld = v;}
  inline void SetTracker(G4bool v=true) {fIsTracker = v;}
  inline void SetCalorimeter(G4bool v=true) {fIsCalorimeter = v;}
  inline G4bool IsCalorimeter() const {return fIsCalorimeter;}

private:
  G4bool fIsWorld;
  G4bool fIsTracker;
  G4bool fIsCalorimeter;
};

#endif
