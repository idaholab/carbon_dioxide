#pragma once

#include "MooseObjectUnitTest.h"
#include "CarbonDioxideLiquidFluidProperties.h"

class CarbonDioxideLiquidFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  CarbonDioxideLiquidFluidPropertiesTest() : MooseObjectUnitTest("CarbonDioxideApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("CarbonDioxideLiquidFluidProperties");
    _fe_problem->addUserObject("CarbonDioxideLiquidFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObjectTempl<CarbonDioxideLiquidFluidProperties>("fp");
  }

  const CarbonDioxideLiquidFluidProperties * _fp;
};
