#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "PhysicsList.hh"

#include "RootIO.hh"


// ----------------------------------------------------------------------
int main(int argc,char** argv) {

  // -- Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // -- Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // -- Default run manager
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif

  // -- Detector construction
  DetectorConstruction* detector = new DetectorConstruction;
  runManager->SetUserInitialization(detector);

  // G4VModularPhysicsList* physicsList = new FTFP_BERT;
  // physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  // runManager->SetUserInitialization(physicsList);
  //  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(new PhysicsList);

  // -- Set user action classes: define particle generator
  runManager->SetUserAction(new PrimaryGeneratorAction(detector));
  runManager->SetUserAction(new RunAction);
  runManager->SetUserAction(new EventAction);

  runManager->Initialize();

  // -- Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  // --  Pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // -- Process macro or start UI session
  if (!ui) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  RootIO::GetInstance()->Close();

  delete visManager;
  delete runManager;

  return 0;
}
