#pragma once

#include "MooseApp.h"

class Factory;
class CarbonDioxideApp;

template <>
InputParameters validParams<CarbonDioxideApp>();

class CarbonDioxideApp : public MooseApp
{
public:
  CarbonDioxideApp(InputParameters parameters);

  static void registerApps();
  static void registerAll(Factory & f, ActionFactory & af, Syntax & s);

protected:
};
