//* This file is part of carbon_dioxide
//* https://github.com/idaholab/carbon_dioxide
//*
//* All rights reserved, see NOTICE.txt for full restrictions
//* https://github.com/idaholab/carbon_dioxide/blob/master/NOTICE.txt
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

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
    _fp = &_fe_problem->getUserObject<CarbonDioxideTwoPhaseFluidProperties>("fp");
  }

  const CarbonDioxideTwoPhaseFluidProperties * _fp;
};
