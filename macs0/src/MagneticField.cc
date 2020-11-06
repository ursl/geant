#include "MagneticField.hh"
#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"

// ----------------------------------------------------------------------
MagneticField::MagneticField() : G4MagneticField(), fMessenger(nullptr), fBy(-0.014*tesla) {
  DefineCommands();
}


// ----------------------------------------------------------------------
MagneticField::~MagneticField() {
  delete fMessenger;
}

// ----------------------------------------------------------------------
void MagneticField::SetField(G4double val) {
  G4cout << "=========> MagneticField::SetField(" << val << ")!" << G4endl;
  fBy = val;
}



// ----------------------------------------------------------------------
void MagneticField::GetFieldValue(const G4double [4], double *bField) const {
  bField[0] = 0.;
  bField[1] = fBy;
  bField[2] = 0.;
  //  G4cout << "=========> MagneticField::GetFieldValue() called, fBy = " << fBy  << G4endl;
}


// ----------------------------------------------------------------------
void MagneticField::DefineCommands() {
  // Define /B5/field command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, "/macs0/field/", "Field control");

  // fieldValue command
  auto& valueCmd = fMessenger->DeclareMethodWithUnit("value","tesla", &MagneticField::SetField, "Set field strength.");
  valueCmd.SetParameterName("field", true);
  valueCmd.SetDefaultValue("-0.014");
}
