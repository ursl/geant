#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(macs0)

#set(CMAKE_CXX_FLAGS "-O0 -g -std=c++17 -fPIC -Wno-return-type")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_BUILD_TYPE Debug)

set(PROJECT_RUN_DIR ${PROJECT_SOURCE_DIR}/run)


#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)

# ROOT version 6 required
if(ROOT_FOUND)
  STRING(REGEX MATCH "6.*" VERSION6MATCH ${ROOT_VERSION})
  if(NOT VERSION6MATCH)
    message(FATAL_ERROR "P01 requires ROOT 6")
  endif()
endif()

# Include ROOT's CMake functions for dictionary generation
#include("${ROOT_DIR}/modules/RootMacros.cmake")


#----------------------------------------------------------------------
FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF()

#----------------------------------------------------------------------------
# P01 requires shared libraries
#
if(NOT Geant4_shared_FOUND)
  message(FATAL_ERROR "P01 must use shared libraries")
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${PROJECT_SOURCE_DIR}/../common/
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS}
		   )
#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc ${PROJECT_SOURCE_DIR}/../common/*.cc)
set(sources
# -- first files from local dir
${PROJECT_SOURCE_DIR}/src/ChamberParameterisation.cc
${PROJECT_SOURCE_DIR}/src/DetectorConstruction.cc
${PROJECT_SOURCE_DIR}/src/DetectorMessenger.cc
${PROJECT_SOURCE_DIR}/src/ElectricFieldSetup.cc
${PROJECT_SOURCE_DIR}/src/EventAction.cc
${PROJECT_SOURCE_DIR}/src/F02ElectricFieldSetup.cc
${PROJECT_SOURCE_DIR}/src/F02FieldMessenger.cc
${PROJECT_SOURCE_DIR}/src/MCPHit.cc
${PROJECT_SOURCE_DIR}/src/MCPSD.cc
${PROJECT_SOURCE_DIR}/src/MagneticField.cc
${PROJECT_SOURCE_DIR}/src/PhysicsList.cc
${PROJECT_SOURCE_DIR}/src/PrimaryGeneratorAction.cc
${PROJECT_SOURCE_DIR}/src/RegionInformation.cc
${PROJECT_SOURCE_DIR}/src/RootIO.cc
${PROJECT_SOURCE_DIR}/src/RunAction.cc
${PROJECT_SOURCE_DIR}/src/StackingAction.cc
${PROJECT_SOURCE_DIR}/src/SteppingAction.cc
${PROJECT_SOURCE_DIR}/src/SteppingVerbose.cc
${PROJECT_SOURCE_DIR}/src/TGenCand.cc
${PROJECT_SOURCE_DIR}/src/TGenVtx.cc
${PROJECT_SOURCE_DIR}/src/THit.cc
${PROJECT_SOURCE_DIR}/src/TrackInformation.cc
${PROJECT_SOURCE_DIR}/src/TrackerHit.cc
${PROJECT_SOURCE_DIR}/src/TrackerSD.cc
${PROJECT_SOURCE_DIR}/src/TrackingAction.cc
${PROJECT_SOURCE_DIR}/src/Trajectory.cc
${PROJECT_SOURCE_DIR}/src/rEvent.cc
# -- second files from common dir
${PROJECT_SOURCE_DIR}/../common/Mu3eMuonInternalConversionDecayWithSpin.cc
${PROJECT_SOURCE_DIR}/../common/Mu3eMuonRadiativeDecayChannelWithSpin.cc
${PROJECT_SOURCE_DIR}/../common/MuDecayChannel.cc
${PROJECT_SOURCE_DIR}/../common/musrMuEnergyLossLandau.cc
${PROJECT_SOURCE_DIR}/../common/musrMuFormation.cc
${PROJECT_SOURCE_DIR}/../common/musrMuStop.cc
${PROJECT_SOURCE_DIR}/../common/musrMuonium.cc
${PROJECT_SOURCE_DIR}/../common/yields.cc
)

#file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)



#----------------------------------------------------------------------------
# Generate dictionaries, add ROOT libraries properties
#
REFLEX_GENERATE_DICTIONARY(Classes include/Classes.hh SELECTION xml/selection.xml)
add_library(ClassesDict SHARED ${sources} Classes.cxx)
set(libsuffix .so)
set(ROOT_LIBRARY_PROPERTIES ${ROOT_LIBRARY_PROPERTIES} SUFFIX ${libsuffix})
set_target_properties(ClassesDict PROPERTIES ${ROOT_LIBRARY_PROPERTIES})
target_link_libraries(ClassesDict  ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(runG4 runG4.cc ${sources})
target_link_libraries(runG4 ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
target_link_libraries(ClassesDict ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(DEPENDS runG4 readHits )

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS runG4 DESTINATION bin)
install(TARGETS ClassesDict DESTINATION lib)
install(FILES ${PROJECT_BINARY_DIR}/Classes_rdict.pcm DESTINATION lib)
