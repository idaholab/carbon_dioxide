//* This file is part of carbon_dioxide
//* https://github.com/idaholab/carbon_dioxide
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/carbon_dioxide/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CarbonDioxideApp.h"
#include "CarbonDioxideRevision.h"
#include "MooseSyntax.h"
#include "AppFactory.h"

// Modules
#include "ModulesApp.h"

InputParameters
CarbonDioxideApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_output_syntax") = false;
  return params;
}

registerKnownLabel("CarbonDioxideApp");

CarbonDioxideApp::CarbonDioxideApp(InputParameters parameters) : MooseApp(parameters)
{
  CarbonDioxideApp::registerAll(_factory, _action_factory, _syntax);
}

// External entry point for dynamic application loading
extern "C" void
CarbonDioxideApp__registerApps()
{
  CarbonDioxideApp::registerApps();
}

void
CarbonDioxideApp::registerApps()
{
  registerApp(CarbonDioxideApp);
}

// External entry point for dynamic object registration
extern "C" void
CarbonDioxideApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  CarbonDioxideApp::registerAll(f, af, s);
}

void
CarbonDioxideApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  Registry::registerObjectsTo(f, {"CarbonDioxideApp"});
  Registry::registerActionsTo(af, {"CarbonDioxideApp"});

  ModulesApp::registerAllObjects<CarbonDioxideApp>(f, af, s);
}
