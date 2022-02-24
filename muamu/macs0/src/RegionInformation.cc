#include "RegionInformation.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RegionInformation::RegionInformation()
  :G4VUserRegionInformation(),
   fIsWorld(false),fIsTracker(false),fIsCalorimeter(false)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RegionInformation::~RegionInformation()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RegionInformation::Print() const
{
 G4cout << "I'm ";
 if(fIsWorld) { G4cout << "World."; }
 else if(fIsTracker) { G4cout << "Tracker."; }
 else if(fIsCalorimeter) { G4cout << "Calorimeter."; }
 else { G4cout << "unknown."; }
 G4cout << G4endl;
}
