//* This file is part of carbon_dioxide
//* https://github.com/idaholab/carbon_dioxide
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/carbon_dioxide/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CarbonDioxideApp.h"
// Moose Includes
#include "MooseInit.h"
#include "Moose.h"
#include "MooseApp.h"
#include "AppFactory.h"

// Create a performance log
PerfLog Moose::perf_log("CarbonDioxide");

// Begin the main program.
int
main(int argc, char * argv[])
{
  // Initialize MPI, solvers and MOOSE
  MooseInit init(argc, argv);

  // Register this application's MooseApp and any it depends on
  CarbonDioxideApp::registerApps();

  std::shared_ptr<MooseApp> app = AppFactory::createAppShared("CarbonDioxideApp", argc, argv);

  // Execute the application
  app->run();

  return 0;
}
