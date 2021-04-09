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
  if (CheckCondition(aStep)) {
    GetDatas(&aStep);
    G4Step theStep;
    if ("Muonium" == particle->GetParticleName()) {
      if (0) G4cout << " -> musrMuFormation::PostStepDoIt> rnd = " << rnd << " -> Muonium formed!!!!!!!!"
		    << G4endl;
    }
    PrepareSecondary(trackData);
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

    fParticleChange.AddSecondary(aSecondary);
    // -- the (original) muon is killed, not the Muonium!
    fParticleChange.ProposeTrackStatus(fStopAndKill) ;
  } else {
    fParticleChange.ProposeTrackStatus(trackData.GetTrackStatus()) ;
  }
  return &fParticleChange;
}


// ----------------------------------------------------------------------
G4bool musrMuFormation::CheckCondition(const G4Step& aStep) {
  G4bool condition=false;
  p_name = aStep.GetTrack()->GetDefinition()->GetParticleName(); // particle name
  std::string logVolName = aStep.GetTrack()->GetVolume()->GetLogicalVolume()->GetName();
  if ((p_name == "mu+") && (logVolName=="Target")) {
      condition=true;
  }
  if (0) G4cout << "condition = " << condition
		<< " ID = " << aStep.GetTrack()->GetTrackID()
		<< " name = " << p_name
		<< " Ekin = " << aStep.GetTrack()->GetKineticEnergy()
		<< G4endl;
  return condition;
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

// ----------------------------------------------------------------------
void musrMuFormation::GetDatas(const G4Step* aStep) {
  particleTable=G4ParticleTable::GetParticleTable();
  rnd=G4UniformRand();
  G4double E = aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy()/CLHEP::keV;
  Gonin.GetYields(E, 105.658369*1000, yvector); // Energy [keV], muon mass [keV/c2], yield table
  G4String p_new = "Muonium";

  if (0) G4cout << "musrMuFormation::GetDatas> rnd = " << rnd
		<< " yvector[0] = " << yvector[0] << " name = " << p_name
		<< G4endl;

  // -- muonium only formed for positive muons entering Target
  if (p_name=="mu+") {
    if (rnd < yvector[0]) {
      particle = particleTable->FindParticle(p_name) ;
      if (0) G4cout << "musrMuFormation::GetDatas> rnd = " << rnd  << " -> stay with muon"
		    << " yvector[0] = " << yvector[0] << " name = " << p_name
		    << G4endl;
    } else {
      particle = particleTable->FindParticle(p_new);
      if (0) G4cout << "musrMuFormation::GetDatas> rnd = " << rnd << " -> Muonium formed!!!!"
		    << " E = " << E
		    << G4endl
		    << " DP E = " << aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy()
		    << "   vec{|p|} = " << aStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection()
		    << G4endl;
    }

    // Set the new dynamic particle DP
    DP = new G4DynamicParticle(particle,
			       aStep->GetTrack()->GetDynamicParticle()->GetMomentumDirection(),
			       aStep->GetTrack()->GetDynamicParticle()->GetKineticEnergy());
    DP->SetProperTime(  aStep->GetTrack()->GetDynamicParticle()->GetProperTime());
    DP->SetPolarization(aStep->GetTrack()->GetDynamicParticle()->GetPolarization().x(),
			aStep->GetTrack()->GetDynamicParticle()->GetPolarization().y(),
			aStep->GetTrack()->GetDynamicParticle()->GetPolarization().z());
    DP->SetPreAssignedDecayProperTime(aStep->GetTrack()->GetDynamicParticle()->GetPreAssignedDecayProperTime());
  }
}

// ----------------------------------------------------------------------
void musrMuFormation::PrepareSecondary(const G4Track& track) {
  if (p_name=="mu+") {
    if (0) cout << "musrMuFormation::PrepareSecondary: name = "
		<< track.GetParticleDefinition()->GetParticleName()
		<< " at x = " << track.GetPosition()
		<< " with p = " << track.GetMomentum() << endl;
    aSecondary = new G4Track(DP,track.GetGlobalTime(),track.GetPosition());
  }
}
