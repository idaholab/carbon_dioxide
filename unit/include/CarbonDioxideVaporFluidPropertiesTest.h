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
