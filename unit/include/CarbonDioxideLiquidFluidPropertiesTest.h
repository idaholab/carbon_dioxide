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
    _fp = &_fe_problem->getUserObject<CarbonDioxideLiquidFluidProperties>("fp");
  }

  const CarbonDioxideLiquidFluidProperties * _fp;
};
