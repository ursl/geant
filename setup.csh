# -- switch to QT based build
setenv GEANT geant4.10.07.p02
setenv GEANT geant4.10.04.p02-withQt-install

# -- subdirectory for examples:
setenv GEANTC Geant4-10.7.2
setenv GEANTC Geant4-10.4.2


cd /Users/ursl/fsx/geant4/${GEANT}/bin && source geant4.csh && cd -
setenv G4INC /Users/ursl/fsx/geant4/${GEANT}/include/Geant4/
setenv G4SRC /Users/ursl/fsx/geant4/${GEANT}-source/source/
setenv G4EX /Users/ursl/fsx/geant4/${GEANT}/share/${GEANTC}/examples/
