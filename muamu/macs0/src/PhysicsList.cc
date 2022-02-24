#include "PhysicsList.hh"

#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"

#include <G4HadronElasticPhysics.hh>
#include <G4HadronPhysicsQGSP_BERT.hh>

#include <G4DecayTable.hh>
#include <G4PhaseSpaceDecayChannel.hh>
#include <G4MuonDecayChannelWithSpin.hh>

#include <G4BaryonConstructor.hh>
#include <G4BosonConstructor.hh>
#include <G4IonConstructor.hh>
#include <G4LeptonConstructor.hh>
#include <G4MesonConstructor.hh>

// -- scattering
#include <G4CoulombScattering.hh>
#include <G4eMultipleScattering.hh>
#include <G4GoudsmitSaundersonMscModel.hh>
#include <G4WentzelVIModel.hh>

// -- process
#include <boost/algorithm/string.hpp>
#include "G4ProcessManager.hh"

#include <G4PhotoElectricEffect.hh>
#include <G4ComptonScattering.hh>
#include <G4GammaConversion.hh>

#include <G4eIonisation.hh>
#include <G4eBremsstrahlung.hh>
#include <G4eplusAnnihilation.hh>

#include <G4MuMultipleScattering.hh>
#include <G4MuIonisation.hh>
#include <G4MuBremsstrahlung.hh>
#include <G4MuPairProduction.hh>

#include <G4MuonVDNuclearModel.hh>
#include <G4MuonNuclearProcess.hh>

#include <G4hMultipleScattering.hh>
#include <G4hIonisation.hh>
#include <G4hBremsstrahlung.hh>
#include <G4hPairProduction.hh>

#include <G4UserSpecialCuts.hh>
#include <G4StepLimiter.hh>

#include <G4DecayWithSpin.hh>

#include <G4BertiniElectroNuclearBuilder.hh>

#include "../common/Mu3eMuonInternalConversionDecayWithSpin.hh"
#include "../common/Mu3eMuonRadiativeDecayChannelWithSpin.hh"

#include "../common/musrMuonium.hh"
#include "../common/MuDecayChannel.hh"
#include "../common/musrMuEnergyLossLandau.hh"
#include "../common/musrMuFormation.hh"
#include "../common/musrMuStop.hh"

// ----------------------------------------------------------------------
PhysicsList::PhysicsList() : G4VUserPhysicsList(){
  SetVerboseLevel(1);

  fPhysicsConstructors.push_back(new G4HadronElasticPhysics(0));
  fPhysicsConstructors.push_back(new G4HadronPhysicsQGSP_BERT(0));

  fEGammaCut = 10.;
  fMuonPolarization = 1.;
}

// ----------------------------------------------------------------------
PhysicsList::~PhysicsList() { }

// ----------------------------------------------------------------------
void PhysicsList::SetCuts() {
  G4VUserPhysicsList::SetCuts();
}

// ----------------------------------------------------------------------
void PhysicsList::ConstructParticle() {
  G4cout << "PhysicsList::ConstructParticle()" << G4endl;

  for(auto pc : fPhysicsConstructors) pc->ConstructParticle();

  G4BaryonConstructor::ConstructParticle();
  G4BosonConstructor::ConstructParticle();
  G4IonConstructor::ConstructParticle();
  G4LeptonConstructor::ConstructParticle();
  G4MesonConstructor::ConstructParticle();

  // -- set up Muonium
  //  musrMuonium::MuoniumDefinition();
  // MuDecayChannel *MuDecay = new MuDecayChannel("Muonium", 1.);

  // -- set up muon decays
  G4double radbr10MeV = 0.014;
  Mu3eMuonRadiativeDecayChannelWithSpin *raddecay
    = new Mu3eMuonRadiativeDecayChannelWithSpin("mu+", radbr10MeV, fEGammaCut);

  G4double radbr = radbr10MeV*raddecay->GetBRScale();
  G4double convbr = 3.4e-5;

  G4VDecayChannel *michel
    = new G4MuonDecayChannelWithSpin("mu+", 1.0-radbr-convbr);

  G4bool isICWeighted(false);

  G4VDecayChannel *intconv
    // = new Mu3eMuonInternalConversionDecayWithSpin("intconv", "mu+", convbr,
    // 						  0, 0 , 0, 0, fMuonPolarization);
    = new Mu3eMuonInternalConversionDecayWithSpin("intconv", "mu+",
						  convbr, 0,
						  1.e7, 0, 0, 0,
						  fMuonPolarization, isICWeighted);

  auto table = new G4DecayTable();
  //  table->Insert(MuDecay); // this is not required?!
  table->Insert(michel);
  table->Insert(raddecay);
  table->Insert(intconv);
  //  table->Insert(MuDecay);
  G4MuonPlus::Definition()->SetDecayTable(table);


}


// ----------------------------------------------------------------------
void PhysicsList::ConstructProcess() {
  G4cout << "PhysicsList::ConstructProcess()" << G4endl;
  AddTransportation();

  for (auto ctor : fPhysicsConstructors) ctor->ConstructProcess();

  struct BertiniElectroNuclearBuilder : G4BertiniElectroNuclearBuilder {
    BertiniElectroNuclearBuilder() : G4BertiniElectroNuclearBuilder(true) {}
    void Build() override {
      G4BertiniElectroNuclearBuilder::Build();
      theModel->SetMaxEnergy(1000 * CLHEP::TeV);
    }
  };

  if(!bertiniElectroNuclearBuilder) {
    bertiniElectroNuclearBuilder = new BertiniElectroNuclearBuilder;
    bertiniElectroNuclearBuilder->Build();
  }

  auto helper = G4PhysicsListHelper::GetPhysicsListHelper();

  auto muonDecay = new G4DecayWithSpin();
  //  muonDecay->SetExtDecayer(&overlapDecayer);

  auto stepLimiter = new G4StepLimiter();

  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ) {
    G4ParticleDefinition* particle = particleIterator->value();
    if (particle->GetProcessManager() == nullptr) {
      G4cout << "PhysicsList::ConstructProcess(): particle->GetParticleName() = "
	     << particle->GetParticleName()
	     << " no process manager" << G4endl;
      continue;
    }
    G4String particleName = particle->GetParticleName();

    if(particleName == "gamma") {
      helper->RegisterProcess(new G4PhotoElectricEffect, particle);
      helper->RegisterProcess(new G4ComptonScattering, particle);
      helper->RegisterProcess(new G4GammaConversion, particle);

      helper->RegisterProcess(new G4UserSpecialCuts, particle);
    }
    else if(particleName == "e+"
	    || particleName == "e-"
	    ) {
      helper->RegisterProcess(new G4CoulombScattering, particle);

      helper->RegisterProcess(new G4eIonisation, particle);
      helper->RegisterProcess(new G4eBremsstrahlung, particle);
      if(particleName == "e+") helper->RegisterProcess(new G4eplusAnnihilation, particle);

      helper->RegisterProcess(new G4UserSpecialCuts, particle);
    }
    else if (particleName == "mu+"
	     || particleName == "mu-"
	     || boost::starts_with(particleName, "mu+/")
	     ) {
      helper->RegisterProcess(new G4CoulombScattering, particle);

      //            helper->RegisterProcess(new G4MuMultipleScattering, particle);
      helper->RegisterProcess(new G4MuIonisation, particle);
      helper->RegisterProcess(new G4MuBremsstrahlung, particle);
      helper->RegisterProcess(new G4MuPairProduction, particle);

      auto muonNuclearProcess = new G4MuonNuclearProcess();
      muonNuclearProcess->RegisterMe(new G4MuonVDNuclearModel);
      helper->RegisterProcess(muonNuclearProcess, particle);

      if (particleName == "mu+") {
	G4ProcessManager* pmanager = particle->GetProcessManager();
	pmanager->AddProcess(new musrMuFormation, -1, -1, 2);
      }
    }
    else if (particleName == "pi+"
	     || particleName == "pi-"
	     || particleName == "proton"
	     ) {
      helper->RegisterProcess(new G4hMultipleScattering, particle);
      helper->RegisterProcess(new G4hIonisation, particle);
      helper->RegisterProcess(new G4hBremsstrahlung, particle);
      helper->RegisterProcess(new G4hPairProduction, particle);
    } else if (particleName == "Muonium") {
      G4ProcessManager* pmanager = particle->GetProcessManager();
      if (1) {
	G4cout << "PhysicsList::ConstructProcess(): particle->GetParticleName() = "
	       << particle->GetParticleName()
	       << " registering  musrMuEnergyLossLandau process "
	       << " with GetVerboseLevel() = " << GetVerboseLevel()
	       << G4endl;

	//helper->RegisterProcess(new musrMuEnergyLossLandau, particle);

	G4VProcess *aMuEnergyLossLandau = new musrMuEnergyLossLandau();
	//	aMuEnergyLossLandau->SetVerboseLevel(GetVerboseLevel());
	pmanager->AddProcess(aMuEnergyLossLandau);
	// -- this is essential for getting this process activated:
	pmanager->SetProcessOrdering(aMuEnergyLossLandau, idxPostStep, 1);
      }
    } else if (!particle->IsShortLived()
	       && particle->GetPDGCharge() != 0
	       && particle->GetParticleName() != "chargedgeantino"
	       ) {
      helper->RegisterProcess(new G4hMultipleScattering, particle);
      helper->RegisterProcess(new G4hIonisation, particle);
    }

    // Add Decay Process
    if(particleName == "mu+"
       || particleName == "mu-"
       || boost::starts_with(particleName, "mu+/")
       ) {
      helper->RegisterProcess(muonDecay, particle);
      G4cout << "PhysicsList::ConstructProcess(): particle->GetParticleName() = "
	     << particle->GetParticleName()
	     << " registering decay processes" << G4endl;

    }
    else {
      auto decay = new G4Decay();
      helper->RegisterProcess(decay, particle);
    }

    // Step limitation seen as a process
    if(particle->GetPDGCharge() != 0) {
      helper->RegisterProcess(stepLimiter, particle);
    }
  }
}
