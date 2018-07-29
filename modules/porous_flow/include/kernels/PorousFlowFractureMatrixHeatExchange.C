//* This file is part of the geothermal

#include "PorousFlowFractureMatrixHeatExchange.h"

#include "MooseVariable.h"

registerMooseObject("PorousFlowApp", PorousFlowFractureMatrixHeatExchange);

template <>
InputParameters
validParams<PorousFlowFractureMatrixHeatExchange>()
{
  InputParameters params = validParams<Kernel>();

  params.addParam<Real>("heat_convection_coef", 0.0, "heat conduction coefficient between fracture fluid and rock surface");
  params.addRequiredCoupledVar(
      "couple_temperature", "The temperature of the matrix or of the fluid will be used.");
  MooseEnum heat_change_type("MatrixToFluid FluidToMatrix", "MatrixToFluid");
  params.addParam<MooseEnum>("heat_change_type", heat_change_type, "the heat flow from matrix to fluid or from fluid to matrxi");
  return params;
}

PorousFlowFractureMatrixHeatExchange::PorousFlowFractureMatrixHeatExchange(const InputParameters & parameters)
  : Kernel(parameters), 
    _heat_change_type(getParam<MooseEnum>("heat_change_type").getEnum<HeatChangeType>()),
    _h(getParam<Real>("heat_convection_coef")),
    _couple_temperature(coupledValue("couple_temperature"))
{
}

Real
PorousFlowFractureMatrixHeatExchange::computeQpResidual()
{
  if (_heat_change_type == HeatChangeType::MatrixToFluid)
    // the couple temperature should be the temperatur of matrix 
    return  _test[_i][_qp] * _h * (_u[_qp] - _couple_temperature[_qp]);
  else if (_heat_change_type == HeatChangeType::FluidToMatrix)
   // the couple temperature should be the temperature of fluid
    return  _test[_i][_qp] * _h * ( _couple_temperature[_qp] - _u[_qp]);
  else
    mooseError("the Heat Change type should be either 'MatrixToFluid' or 'FluidToMatrux'");
}

Real
PorousFlowFractureMatrixHeatExchange::computeQpJacobian()
{
  return _test[_i][_qp]* _h * _phi[_j][_qp];
}
