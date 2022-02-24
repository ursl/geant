/// \file Mu3eMuonInternalConversionDecayWithSpin.h

/*
 * Mu3eInternalConversionDecay.h
 *
 *  Created on: Sept 20, 2012
 *      Author: nberger
 */

#ifndef Mu3eMuonInternalConversionDecayWithSpin_H_
#define Mu3eMuonInternalConversionDecayWithSpin_H_

#include "G4VDecayChannel.hh"

class TH1D;
class TH2D;

/// Internal conversion decay of the muon
/**
 *  Generates a five particle phase space configuration with RAMBO,
 * 	then does accept/reject with the matrix element by
 * 	R. M. Djilkibaev and R. V. Konoplich, “Rare Muon Decay mu+ -> e+e-e+ nu_mu nu_e”,
 *  Phys.Rev., D79 073004, 2009.
 *
 *  Additional cuts on the decay electron kinematics can be specified in the constructor
 *
 *  Currently uses an unpolarized matrix element but is a polarized decay for Geant4 in order
 *  to allow mixing with polarized matrix elements (e.g. for the radiative decay)
 *
 */
struct Mu3eMuonInternalConversionDecayWithSpin : G4VDecayChannel {
  Mu3eMuonInternalConversionDecayWithSpin(
					  const G4String& theName,                ///< Name of the decay channel
					  const G4String& theParentName,          ///< Name of the parent particle
					  G4double        theBR,                  ///< Branching fraction of the decay
					  G4int           Selection,              ///< Type of selection; 0 cut on all three electrons, 1 cut on an e+e- pair
					  G4double        CosThetaCut,            ///< Cut on the absolute value of the polar angle of the decay electrons, 0 for no cut
					  G4double        EminCut,                ///< Cut on the minimum energy of the decay electrons, 0 for no cut
					  G4double        EvisCut,                ///< Cut on the minimum energy visible in the three electrons, 0 for no cut
					  G4double        MassCut,                ///< Cut on the minimum invariant mass in the three electrons, 0 for no cut
					  G4double        theParentPolarization,  ///< Polarization of the parent particle
					  G4bool          weighted = false        ///< Instead of accept/reject, generate event weights
					  );

  /// Perform the actual decay
  virtual G4DecayProducts *DecayIt(G4double);

  /// Set the polarization of the parent particle, which is currently ignored
  void SetPolarization(G4ThreeVector polar){ parent_polarization = polar;};
  /// Get the polarization of the parent particle
  const G4ThreeVector& GetPolarization() const {return parent_polarization;};
  /// Get the branching fraction for the current cuts and fill the histograms
  double GetBranchingFraction(double nevents, TH1D * allmasshisto, TH1D * selmasshisto, TH1D * espectrum, TH1D *eminspectrum, TH2D * eespectrum, TH2D * eespectrumsel, TH2D * dalitz, TH1D * etothisto, TH1D * costhposhisto, TH1D * costhelhhisto, TH1D * costhellhisto, TH1D * costhposhistosel, TH1D * costhelhhistosel, TH1D * costhellhistosel, TH1D * ptposhistosel, TH1D * ptelhhistosel, TH1D * ptellhistosel);


  /// Perform selection on the generated events, return e+e-e+ invariant mass
  int selection(double P[5][4], double & mass);
  /// Perform selection on the generated events
  int selection(double P[5][4]);
  /// Perform selection on the generated events, return e+e-e+ invariant mass and cos theta of the electrons
  int selection(double P[5][4], double & mass, double & costhpos, double & costhelh, double & costhell);
  /// Perform selection on the generated events, return e+e-e+ invariant mass and cos theta of the electrons, cut on pt
  int selection(double P[5][4], double & mass, double & costhpos, double & costhelh, double & costhell, double ptcut, double eviscut, double & ptelh, double & ptell, double & ptpos, double & evis);
  /// Matrix element
  double matmu3e2nu(double P[5][4]);
  double matmu3e2nu_DK(double P[5][4]);       // previous matrix element by Djilkibaev & Konoplich
protected:
  /// Polarization vector of the muon
  G4ThreeVector parent_polarization;

  /// Type of selection; 0 cut on all three electrons, 1 cut on an e+e- pair
  const G4int    selectionType;
  /// Cut on the absolute value of the polar angle of the decay electrons, 0 for no cut
  const G4double fCosThetaCut;
  /// Cut on the minimum energy of the decay electrons, 0 for no cut
  const G4double fEminCut;
  /// Cut on the minimum energy visible in the three electrons, 0 for no cut
  const G4double fEvisCut;
  /// Cut on the minimum invariant mass in the three electrons, 0 for no cut
  const G4double fMassCut;
  /// Instead of accept/reject, generate event weights
  const G4bool weighted;
  /// Maximum matrix element
  G4double matmax;

};

#endif /* MU3EOVERLAPDECAY_H_ */
