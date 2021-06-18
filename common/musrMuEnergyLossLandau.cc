/***************************************************************************
 *  musrSim - the program for the simulation of (mainly) muSR instruments. *
 *          More info on http://lmu.web.psi.ch/simulation/index.html .     *
 *          musrSim is based od Geant4 (http://geant4.web.cern.ch/geant4/) *
 *                                                                         *
 *  Copyright (C) 2009 by Paul Scherrer Institut, 5232 Villigen PSI,       *
 *                                                       Switzerland       *
 *                                                                         *
 *  This program is free software; you can redistribute it and/or modify   *
 *  it under the terms of the GNU General Public License as published by   *
 *  the Free Software Foundation; either version 2 of the License, or      *
 *  (at your option) any later version.                                    *
 *                                                                         *
 *  This program is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *  GNU General Public License for more details.                           *
 *                                                                         *
 *  You should have received a copy of the GNU General Public License      *
 *  along with this program; if not, write to the Free Software            *
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              *
 ***************************************************************************/

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//  Muon energy loss in thin C-foil adding a Landau distribution to the energy
//  loss
//  Id    : musrMuEnergyLossLandau, v 1.0
//  Author: Thomas Prokscha
//  Date  : 2016-08
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#include "musrMuEnergyLossLandau.hh"

using namespace std;

musrMuEnergyLossLandau::musrMuEnergyLossLandau(const G4String& name, G4ProcessType aType)
               : G4VDiscreteProcess(name, aType)
{
  random = new TRandom();
}

musrMuEnergyLossLandau::~musrMuEnergyLossLandau(){}

// double musrMuEnergyLossLandau::landauMPV   = 0.1;
// double musrMuEnergyLossLandau::landauSigma = 0.1;

double musrMuEnergyLossLandau::landauMPV   = 0.04;
double musrMuEnergyLossLandau::landauSigma = 0.02;

// ----------------------------------------------------------------------
G4VParticleChange* musrMuEnergyLossLandau::PostStepDoIt(const G4Track& trackData, const G4Step&  aStep) {
  fParticleChange.Initialize(trackData);
  //  G4cout << "musrMuEnergyLossLandau::PostStepDoIt" << G4endl;

  // -- check whether Mu hit the endplate. If yes, stop and decay it there.
  G4String p_name = trackData.GetDefinition()->GetParticleName(); // particle name
  if (p_name != "Muonium") {
    return &fParticleChange;
  }

  std::string logVolName = trackData.GetVolume()->GetLogicalVolume()->GetName();
  if (logVolName == "Endplate") {
    fParticleChange.Initialize(trackData);
    //    fParticleChange.ProposeTrackStatus(fStopAndKill) ;
    fParticleChange.ProposeTrackStatus(fStopButAlive) ;
    return &fParticleChange;
  }

  // -- unwrapped:  if (CheckCondition(aStep)) {
  if (logVolName == "Target") {

    // -- unwrapped:  GetFinalEnergy(&aStep);
    particleTable = G4ParticleTable::GetParticleTable();
    G4double Eloss, Efinal;
    G4double E = trackData.GetDynamicParticle()->GetKineticEnergy()/CLHEP::keV;

    do {
      Eloss = random->Landau(landauMPV,landauSigma);
      if (Eloss > 0.) break;
    } while (1);

    Efinal = E - Eloss;
    if (Efinal < 0. ) Efinal = 0.;
    if (0) G4cout << "musrMuEnergyLossLandau::GetFinalEnergy: "
		  << " E, Eloss, Efinal = " << E << ", " << Eloss << ", "  << Efinal
		  << G4endl;

    particle = particleTable->FindParticle(p_name) ;

    // Set the new dynamic particle DP
    DP = new G4DynamicParticle(particle,
			       trackData.GetDynamicParticle()->GetMomentumDirection(),
			       Efinal/1000.);

    DP->SetProperTime(  trackData.GetDynamicParticle()->GetProperTime());
    DP->SetPolarization(trackData.GetDynamicParticle()->GetPolarization().x(),
			trackData.GetDynamicParticle()->GetPolarization().y(),
			trackData.GetDynamicParticle()->GetPolarization().z());
    DP->SetPreAssignedDecayProperTime(trackData.GetDynamicParticle()->GetPreAssignedDecayProperTime());

    //G4Step theStep;
    // NOTE: This is where I should change everything!
    //       - model this according to
    // -- unwrapped:  PrepareSecondary(trackData);
    aSecondary = new G4Track(DP, trackData.GetGlobalTime(), trackData.GetPosition());

    fParticleChange.AddSecondary(aSecondary);
    fParticleChange.ProposeTrackStatus(fStopAndKill) ;
  } else {
    fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;
  }
  return &fParticleChange;
}


// ----------------------------------------------------------------------
G4double musrMuEnergyLossLandau::GetMeanFreePath(const G4Track&t,  G4double, G4ForceCondition* condition1) {
  *condition1 = Forced;
  // std::string logVolName = t.GetVolume()->GetLogicalVolume()->GetName();
  // if (logVolName == "Target") {
  //   return 2*CLHEP::nm;
  // }
  return DBL_MAX;
}
