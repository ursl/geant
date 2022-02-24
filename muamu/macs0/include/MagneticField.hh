#ifndef MagneticField_H
#define MagneticField_H

#include "G4MagneticField.hh"

class G4GenericMessenger;

// ----------------------------------------------------------------------
class MagneticField: public G4MagneticField {
public:

  MagneticField();
  ~MagneticField();

  virtual void GetFieldValue(const G4double point[4],double* bField ) const;
  void SetField(G4double val);
  G4double GetField() const { return fBy; }

private:
  void DefineCommands();

  G4GenericMessenger* fMessenger;
  G4double fBy;

};

#endif
