//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file main.cc
/// \brief Main program of the coincidence event

#include "iostream"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
#include "random"
#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "PhysicsList.hh"
#include "time.h"
#include "yaml-cpp/yaml.h"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


int main(int argc,char** argv)
{
  char* fileConfig= argv[1];
  YAML::Node rootNode = YAML::LoadFile(fileConfig);
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
  //G4long seed = time(NULL);//
  //G4long seed = 1128;
  std::random_device se;
  CLHEP::HepRandom::setTheSeed(se());
  //CLHEP::HepRandom::setTheSeed(1128);
  if (argc != 3) {
    std::cerr << "Wrong input arguments! The right way is: ./main [config/config.yaml] [config/optical_property.yaml]";
    return -1;
  }
  
  G4RunManager* runManager = new G4RunManager;
  // Set mandatory initialization classes
  //
  runManager->SetUserInitialization(new DetectorConstruction(argv[1],argv[2]));

  runManager->SetUserInitialization(new PhysicsList());
    
  // Set user action classes
  runManager->SetUserInitialization(new ActionInitialization());
  
  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/tracking/verbose 0");
  UImanager->ApplyCommand("/run/verbose      0");
  UImanager->ApplyCommand("/control/verbose  0");
  UImanager->ApplyCommand("/analysis/verbose 0");
  UImanager->ApplyCommand("/event/verbose    0");
  UImanager->ApplyCommand("/process/verbose  0");
  runManager->Initialize(); 
  //runManager->BeamOn(std::stoi(argv[3]));
  G4long beamon =rootNode["BeamOn"].as<long>();
  runManager->BeamOn(beamon);
  delete runManager;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
