#ifndef   musrMuFormation_h
#define   musrMuFormation_h 1

#include "G4VDiscreteProcess.hh"
#include "G4ParticleTable.hh"

#include "yields.hh"

/*! musrMuFormation class defines the muonium formation process in the Carbon foil
 *  according to yields from Gonin's paper Sci. Rev. Instrum. 65(3), 648-652 (1994).
 * \image html yields3.gif The muonium formation yields.
 * The main parameters are the foil thickness and muon energy. For a given energy,
 * a corresponding proportion of the muons will be converted into Muonium.
 * Concretely, the muon is eliminated and replaced by a Muonium with identical
 * properties, including time, energy, momentum, position etc.
 *
 * The process is executed at the END of a step, i.e. the muon is converted into
 * Muonium AFTER flying through the Carbon foil (see also yields.hh). */

class musrMuFormation : public G4VDiscreteProcess {
 public:

   musrMuFormation(const G4String& name = "musrMuFormation", // process description
		   G4ProcessType aType = fElectromagnetic);

  ~musrMuFormation();


  G4bool IsApplicable(const G4ParticleDefinition& particle);

  //! - Main method. Muonium formation process is executed at the END of a step. */
  G4VParticleChange* PostStepDoIt(
			     const G4Track&,
			     const G4Step&);

  G4double GetMeanFreePath(const G4Track& aTrack,
			   G4double previousStepSize,
			   G4ForceCondition* condition);


  Yields            fGonin;
  G4VParticleChange fParticleChange;
};

#endif
