#pragma once

#include "MooseObjectUnitTest.h"
#include "CarbonDioxideHEMFluidProperties.h"

class CarbonDioxideHEMFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  CarbonDioxideHEMFluidPropertiesTest() : MooseObjectUnitTest("CarbonDioxideApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("CarbonDioxideHEMFluidProperties");
    _fe_problem->addUserObject("CarbonDioxideHEMFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<CarbonDioxideHEMFluidProperties>("fp");
  }

  const CarbonDioxideHEMFluidProperties * _fp;
};
