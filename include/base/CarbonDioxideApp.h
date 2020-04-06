#pragma once

#include "MooseApp.h"

class Factory;

class CarbonDioxideApp : public MooseApp
{
public:
  CarbonDioxideApp(InputParameters parameters);

public:
  static InputParameters validParams();
  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);
};
