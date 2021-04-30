#ifndef Trajectory_h
#define Trajectory_h 1

#include "G4VTrajectory.hh"
#include "G4Allocator.hh"
#include <stdlib.h>
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4TrajectoryPoint.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include <vector>

class G4Polyline;
class G4AttDef;
class G4AttValue;

typedef std::vector<G4VTrajectoryPoint*> TrajectoryPointContainer;

class Trajectory : public G4VTrajectory {
public:
  Trajectory(const G4Track* aTrack);
  virtual ~Trajectory();

  virtual void ShowTrajectory(std::ostream& os=G4cout) const;
  virtual void DrawTrajectory() const;
  virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
  virtual std::vector<G4AttValue>* CreateAttValues() const;
  virtual void AppendStep(const G4Step* aStep);
  virtual void MergeTrajectory(G4VTrajectory* secondTrajectory);

  inline void* operator new(size_t);
  inline void  operator delete(void*);
  inline int operator == (const Trajectory& right) const  {return (this==&right);}

  virtual G4int GetTrackID() const { return fTrackID; }
  virtual G4int GetParentID() const { return fParentID; }
  virtual G4String GetParticleName() const { return fParticleName; }
  virtual G4double GetCharge() const { return fPDGCharge; }
  virtual G4double GetMass() const { return fPDGMass; }
  virtual G4int GetPDGEncoding() const { return fPDGEncoding; }
  virtual G4ThreeVector GetInitialMomentum() const { return fMomentum; }
  virtual int GetPointEntries() const { return fPositionRecord->size(); }
  virtual G4VTrajectoryPoint* GetPoint(G4int i) const { return (*fPositionRecord)[i]; }
  virtual G4double GetLocalTime() const { return fLocalTime; }
  virtual G4double GetGlobalTime() const { return fGlobalTime; }
  virtual G4double GetKineticEnergy() const { return fKineticEnergy; }

private:
  TrajectoryPointContainer* fPositionRecord;
  G4int                        fTrackID;
  G4int                        fParentID;
  G4int                        fTrackStatus;
  G4ParticleDefinition*        fParticleDefinition;
  G4String                     fParticleName;
  G4double                     fPDGCharge;
  G4double                     fPDGMass;
  G4int                        fPDGEncoding;
  G4ThreeVector                fMomentum;
  G4double                     fKineticEnergy;
  G4ThreeVector                fVertexPosition;
  G4double                     fGlobalTime;
  G4double                     fLocalTime;

};

extern G4ThreadLocal G4Allocator<Trajectory> * myTrajectoryAllocator;

inline void* Trajectory::operator new(size_t) {
  if(!myTrajectoryAllocator)
    myTrajectoryAllocator = new G4Allocator<Trajectory>;
  return (void*)myTrajectoryAllocator->MallocSingle();
}

inline void Trajectory::operator delete(void* aTrajectory) {
  myTrajectoryAllocator->FreeSingle((Trajectory*)aTrajectory);
}

#endif
