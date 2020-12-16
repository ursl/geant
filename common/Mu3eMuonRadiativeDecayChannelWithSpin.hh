//
// ------------------------------------------------------------
//      Mu3eMuonRadiativeDecayChannelWithSpin.hh
//
//	Created on :	Oct 10, 2014
//	Author :	aperrevoort
//
//               Samples Radiative Muon Decay
//               References:
//                    TRIUMF/TWIST Technote TN-55:
//                    "Radiative muon decay" by P. Depommier and A. Vacheret
//                    ------------------------------------------------------
//                    Yoshitaka Kuno and Yasuhiro Okada
//                    "Muon Decays and Physics Beyond the Standard Model"
//                    Rev. Mod. Phys. 73, 151 (2001)
//
// ------------------------------------------------------------
//
//
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef MU3EMUONRADIATIVEDECAYCHANNELWITHSPIN_H_
#define MU3EMUONRADIATIVEDECAYCHANNELWITHSPIN_H_

#include "G4VDecayChannel.hh"

struct Mu3eMuonRadiativeDecayChannelWithSpin : G4VDecayChannel {
  // Class Decription
  // This class describes radiative muon decay kinemtics, but
  // gives incorrect energy spectrum for neutrinos

    Mu3eMuonRadiativeDecayChannelWithSpin(
        const G4String& theParentName,
        double          the10MeVBR,
        double          theRadGammaCut
    );

    virtual G4DecayProducts *DecayIt(double theRadGammaCut);

    void SetPolarization(const G4ThreeVector& p) { parent_polarization = p; }
    const G4ThreeVector& GetPolarization() const { return parent_polarization; }

    double GetBRScale() {return scale;}

protected:
    double ScaleBR(double gammacut);

    G4ThreeVector parent_polarization;

    double fron(double Pmu, double x, double y,
                double cthetaE, double cthetaG, double cthetaEG);

    double kuno(double Pmu, double x, double y,
                double cthetaE, double cthetaG, double cthetaEG);

    double scale;
    double gammacut;

};

#endif /* MU3ERADIATIVEDECAYCHANNELWITHSPIN_H_ */
