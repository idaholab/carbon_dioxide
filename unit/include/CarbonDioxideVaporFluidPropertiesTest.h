#pragma once

#include "MooseObjectUnitTest.h"
#include "CarbonDioxideVaporFluidProperties.h"

class CarbonDioxideVaporFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  CarbonDioxideVaporFluidPropertiesTest() : MooseObjectUnitTest("CarbonDioxideApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("CarbonDioxideVaporFluidProperties");
    _fe_problem->addUserObject("CarbonDioxideVaporFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<CarbonDioxideVaporFluidProperties>("fp");
  }

  const CarbonDioxideVaporFluidProperties * _fp;
};
