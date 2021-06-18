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
//  Muonium Formation according to yield.cc function (through GetYields method).
//  Id    : musrMuFormation.cc, v 1.4
//  Author: Taofiq PARAISO, T. Shiroka
//  Date  : 2007-12
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#include "musrMuFormation.hh"

#include "RootIO.hh"

using namespace std;

// ----------------------------------------------------------------------
musrMuFormation::musrMuFormation(const G4String& name, G4ProcessType aType)
               : G4VDiscreteProcess(name, aType){}

// ----------------------------------------------------------------------
musrMuFormation::~musrMuFormation(){}

// ----------------------------------------------------------------------
G4VParticleChange* musrMuFormation::PostStepDoIt(const G4Track& trackData, const G4Step&  aStep) {
  G4VParticleChange fParticleChange;
  fParticleChange.Initialize(trackData);

  G4Track theNewTrack;
  if (0)   G4cout << "START "
		  << " ID = " << trackData.GetTrackID()
		  << " " << trackData.GetDefinition()->GetParticleName()
		  << " Ekin [keV] = " << trackData.GetKineticEnergy()/CLHEP::keV
		  << " L [mm] = " << trackData.GetTrackLength()/CLHEP::mm
		  << " Edep[keV] = " << aStep.GetTotalEnergyDeposit()/CLHEP::keV
		  << " in " << aStep.GetTrack()->GetVolume()->GetLogicalVolume()->GetName()
		  << G4endl;

  // -- muonium only formed for positive muons entering Target
  G4String  p_name = trackData.GetDefinition()->GetParticleName(); // particle name
  if (p_name != "mu+") {
    fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;
    return &fParticleChange;
  }
  std::string logVolName = trackData.GetVolume()->GetLogicalVolume()->GetName();
  if (logVolName != "Target") {
    fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;
    return &fParticleChange;
  }


  // --GetDatas(&aStep);
  G4ParticleDefinition *particle;
  G4DynamicParticle *DP;
  G4double yvector[3];
  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  G4double rnd = G4UniformRand();
  G4double E = trackData.GetDynamicParticle()->GetKineticEnergy()/CLHEP::keV;
  fGonin.GetYields(E, 105.658369*1000, yvector); // Energy [keV], muon mass [keV/c2], yield table

  if (0) G4cout << "musrMuFormation::GetDatas> rnd = " << rnd
		<< " yvector[0] = " << yvector[0] << " name = " << p_name
		<< G4endl;

  if (rnd < yvector[0]) {
    particle = particleTable->FindParticle(p_name) ;
    if (0) G4cout << "musrMuFormation::GetDatas> rnd = " << rnd  << " -> stay with muon"
		  << " yvector[0] = " << yvector[0] << " name = " << p_name
		  << G4endl;
  } else {
    particle = particleTable->FindParticle("Muonium");
    if (0) G4cout << "musrMuFormation::GetDatas> rnd = " << rnd << " -> Muonium formed!!!!"
		  << " E = " << E
		  << G4endl
		  << " DP E = " << trackData.GetDynamicParticle()->GetKineticEnergy()
		  << "   vec{|p|} = " << trackData.GetDynamicParticle()->GetMomentumDirection()
		  << G4endl;
  }
  // Set the new dynamic particle DP
  DP = new G4DynamicParticle(particle,
			     trackData.GetDynamicParticle()->GetMomentumDirection(),
			     trackData.GetDynamicParticle()->GetKineticEnergy());
  DP->SetProperTime(  trackData.GetDynamicParticle()->GetProperTime());
  DP->SetPolarization(trackData.GetDynamicParticle()->GetPolarization().x(),
		      trackData.GetDynamicParticle()->GetPolarization().y(),
		      trackData.GetDynamicParticle()->GetPolarization().z());
  DP->SetPreAssignedDecayProperTime(trackData.GetDynamicParticle()->GetPreAssignedDecayProperTime());
  G4Track *aSecondary = new G4Track(DP, trackData.GetGlobalTime(), trackData.GetPosition());

  // -- PrepareSecondary(trackData);
  if (1) {
    RootIO *rio = RootIO::GetInstance();
    TGenVtx *pVtx = rio->getEvent()->addGenVtx();
    pVtx->fNumber = rio->getEvent()->nGenVtx()-1;
    pVtx->fID = -1313;
    pVtx->fvMom.push_back(trackData.GetTrackID());
    pVtx->fvDau.push_back(aSecondary->GetTrackID());
    pVtx->fGlobalTime = trackData.GetGlobalTime();
    pVtx->fLocalTime = trackData.GetLocalTime();
    pVtx->fV.SetX(trackData.GetPosition().x());
    pVtx->fV.SetY(trackData.GetPosition().y());
    pVtx->fV.SetZ(trackData.GetPosition().z());
  }

  fParticleChange.AddSecondary(aSecondary);
  // -- the (original) muon is killed, not the Muonium!
  fParticleChange.ProposeTrackStatus(fStopAndKill) ;

  return &fParticleChange;
}


// ----------------------------------------------------------------------
G4double musrMuFormation::GetMeanFreePath(const G4Track&t, G4double, G4ForceCondition* condition) {
  *condition = Forced;
  return DBL_MAX;

  std::string logVolName = t.GetVolume()->GetLogicalVolume()->GetName();
  if (logVolName == "Target") {
    return 2*CLHEP::nm;
  }
}
