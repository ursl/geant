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
#include "TrackingAction.hh"
#include "StackingAction.hh"
#include "SteppingAction.hh"
#include "PhysicsList.hh"

#include "RootIO.hh"

using namespace std;

// ----------------------------------------------------------------------
int main(int argc, char** argv) {

  // -- command line arguments
  bool startGUI(false);
  int nevt(100), verbose(0);
  float sgE(-1.), bgE(-1.);
  int sgN(-1.), bgN(-1.);
  G4String cmd("vis.mac"), filename("nada");
  for (int i = 0; i < argc; i++){
    if (!strcmp(argv[i], "-f"))    {filename = argv[++i];}
    if (!strcmp(argv[i], "-g"))    {startGUI = true;}
    if (!strcmp(argv[i], "-c"))    {cmd = argv[++i];}
    if (!strcmp(argv[i], "-n"))    {nevt = atoi(argv[++i]);}
    if (!strcmp(argv[i], "-sgE"))  {sgE  = atof(argv[++i]);}  // energy of sg particles
    if (!strcmp(argv[i], "-bgE"))  {bgE  = atof(argv[++i]);}  // energy of bg particles
    if (!strcmp(argv[i], "-sgN"))  {sgN  = atoi(argv[++i]);}  // number of sg particles
    if (!strcmp(argv[i], "-bgN"))  {bgN  = atoi(argv[++i]);}  // number of bg particles
    if (!strcmp(argv[i], "-v"))    {verbose = atoi(argv[++i]);}
  }

  if ("nada" == filename) {
    filename = cmd;
    size_t pos1 = filename.find(".mac");
    filename = filename.substr(0, pos1);
    filename += ".root";
  }
  cout << "filename: " << filename <<endl;

  RootIO *rio = RootIO::GetInstance(filename);
  if (verbose > 0) rio->setVerbose(verbose);

  // -- Detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = 0;
  if (startGUI) {
    ui = new G4UIExecutive(argc, argv);
  }

  // -- Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // -- Default run manager
#ifdef G4MULTITHREADED
  cout << "runG4> G4MULTITHREADED" << endl;
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  cout << "runG4> G4 single threaded" << endl;
  G4RunManager* runManager = new G4RunManager;
#endif

  // -- Detector construction
  DetectorConstruction* detector = new DetectorConstruction;
  if (verbose > 0) detector->setVerbose(verbose);
  detector->parseCmd(cmd);

  runManager->SetUserInitialization(detector);

  // G4VModularPhysicsList* physicsList = new FTFP_BERT;
  // physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  // runManager->SetUserInitialization(physicsList);
  //  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  PhysicsList *pl = new PhysicsList;
  if (verbose > 0) {
    pl->SetVerboseLevel(verbose);
  }
  runManager->SetUserInitialization(pl);

  // -- Set user action classes: define particle generator
  runManager->SetUserAction(new PrimaryGeneratorAction(detector));
  runManager->SetUserAction(new RunAction);
  runManager->SetUserAction(new EventAction);

  // Test:
  runManager->SetUserAction(new StackingAction);
  runManager->SetUserAction(new SteppingAction);
  runManager->SetUserAction(new TrackingAction);

  if (verbose > 0) {
    ((EventAction*)runManager->GetUserEventAction())->setVerbose(verbose);
    ((RunAction*)runManager->GetUserRunAction())->setVerbose(verbose);
    ((PrimaryGeneratorAction*)runManager->GetUserPrimaryGeneratorAction())->setVerbose(verbose);

    ((StackingAction*)runManager->GetUserStackingAction())->setVerbose(verbose);
  }

  ((EventAction*)runManager->GetUserEventAction())->setNevt(nevt);


  //  runManager->Initialize();

  // --  Pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  G4VisManager* visManager(0);
  // -- Process macro or start (G)UI session
  if (!ui) {
    // batch mode: after *.mac the program exits
    G4String command = "/control/execute ";
    G4String fileName = cmd;
    UImanager->ApplyCommand(command+fileName);
    // -- override default values in macro after executing it
    if (bgE > -0.9) UImanager->ApplyCommand(Form("/macs0/generator/bgKinEnergy %7.6f MeV", bgE));
    if (sgE > -0.9) UImanager->ApplyCommand(Form("/macs0/generator/sgKinEnergy %7.6f MeV", sgE));
    // -- run some events
    UImanager->ApplyCommand(Form("/run/beamOn %d", nevt));
  } else {
    // -- Initialize visualization
    visManager = new G4VisExecutive;
    visManager->Initialize();

    // GUI mode: after vis.mac the program waits for user input
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  }

  RootIO::GetInstance()->Close();

  if (visManager) delete visManager;
  delete runManager;

  return 0;
}
