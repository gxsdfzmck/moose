//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef CO2FLUIDPROPERTIESTEST_H
#define CO2FLUIDPROPERTIESTEST_H

#include "MooseObjectUnitTest.h"
#include "CO2FluidProperties.h"

class CO2FluidPropertiesTest : public MooseObjectUnitTest
{
public:
  CO2FluidPropertiesTest() : MooseObjectUnitTest("MooseUnitApp")
  {
    registerObjects(_factory);
    buildObjects();
  }

protected:
  void registerObjects(Factory & factory) { registerUserObject(CO2FluidProperties); }

  void buildObjects()
  {
    InputParameters uo_pars = _factory.getValidParams("CO2FluidProperties");
    _fe_problem->addUserObject("CO2FluidProperties", "fp", uo_pars);
    _fp = &_fe_problem->getUserObject<CO2FluidProperties>("fp");
  }

  const CO2FluidProperties * _fp;
};

#endif // CO2FLUIDPROPERTIESTEST_H
