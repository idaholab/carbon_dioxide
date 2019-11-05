#pragma once

#include "MooseObjectUnitTest.h"
#include "CarbonDioxideTwoPhaseFluidProperties.h"

class CarbonDioxideTwoPhaseFluidPropertiesTest : public MooseObjectUnitTest
{
public:
  CarbonDioxideTwoPhaseFluidPropertiesTest() : MooseObjectUnitTest("CarbonDioxideApp")
  {
    buildObjects();
  }

protected:
  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("CarbonDioxideTwoPhaseFluidProperties");
    _fe_problem->addUserObject("CarbonDioxideTwoPhaseFluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObjectTempl<CarbonDioxideTwoPhaseFluidProperties>("fp");
  }

  const CarbonDioxideTwoPhaseFluidProperties * _fp;
};
