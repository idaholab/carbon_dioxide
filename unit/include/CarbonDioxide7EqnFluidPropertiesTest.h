#pragma once

#include "MooseObjectUnitTest.h"
#include "CarbonDioxide7EqnFluidProperties.h"

class CarbonDioxide7EqnFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  CarbonDioxide7EqnFluidPropertiesTest() : MooseObjectUnitTest("CarbonDioxideApp") { buildObjects(); }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("CarbonDioxide7EqnFluidProperties");
    _fe_problem->addUserObject("CarbonDioxide7EqnFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObjectTempl<CarbonDioxide7EqnFluidProperties>("fp");
  }

  const CarbonDioxide7EqnFluidProperties * _fp;
};
