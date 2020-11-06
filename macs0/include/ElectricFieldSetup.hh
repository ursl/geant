#ifndef ElectricFieldSetup_h
#define ElectricFieldSetup_h 1

#include "G4ElectricField.hh"
#include "G4UniformElectricField.hh"

class G4FieldManager;
class G4ChordFinder;
class G4EquationOfMotion;
class G4Mag_EqRhs;
class G4EqMagElectricField;
class G4MagIntegratorStepper;
class G4MagInt_Driver;
class G4GenericMessenger;

class ElectricFieldSetup {
public:
  ElectricFieldSetup(G4ThreeVector);  //  The value of the field
  ElectricFieldSetup();               //  A zero field - true value set later

  virtual ~ElectricFieldSetup();

   // Methods to set parameters or select
  void SetStepperType( G4int i) { fStepperType = i ; CreateStepper(); }

  void SetMinStep(G4double s) { fMinStep = s ; }

  void SetFieldValue(G4ThreeVector fieldVector);
  void SetConstantFieldValue(G4double      fieldValue);
  G4ThreeVector GetConstantFieldValue();
   // Set/Get Field strength in Geant4 units

  void UpdateIntegrator();
   // Prepare all the classes required for tracking - from stepper
   //    to Chord-Finder
   //   NOTE:  field and equation must have been created before calling this.

protected:

  // Find the global Field Manager

  G4FieldManager*         GetGlobalFieldManager();

  void CreateStepper();
   // Implementation method - should not be exposed

private:
  G4double                fMinStep;
  G4bool                  fVerbose;

  G4FieldManager*         fFieldManager;

  G4ChordFinder*          fChordFinder;

  G4EqMagElectricField*   fEquation;

  G4ElectricField*        fEMfield;

  G4ThreeVector           fElFieldValue;

  G4MagIntegratorStepper* fStepper;
  G4MagInt_Driver*        fIntgrDriver;

  G4int                   fStepperType;


  void DefineCommands();
  G4GenericMessenger* fMessenger;

};

#endif
