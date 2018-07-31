//*

#include "PorousFlowThermalConductivityTensor.h"

registerMooseObject("PorousFlowApp", PorousFlowThermalConductivityTensor);

template <>
InputParameters
validParams<PorousFlowThermalConductivityTensor>()
{
  InputParameters params = validParams<PorousFlowMaterialVectorBase>();
  params.addRequiredParam<RealTensor>("rock_thermal_conductivity", "The thermal conductivity of rock");
  params.addRequiredParam<RealTensor>("fluid_thermal_conductivity", "The thermal conductivity of fluid");
  params.addParam<Real>("power", 1.5, "the power component of Bruggeman correction");
  params.addParam<bool>("modify", false, "If true, thermal conductivity wiil be modified.");  
  params.addClassDescription("Provide tensor formed thermal conductivity for matrix and fluid");
  return params;
}

PorousFlowThermalConductivityTensor::PorousFlowThermalConductivityTensor(const InputParameters & parameters)
  : PorousFlowMaterialVectorBase(parameters),
    _lambda_rock(getParam<RealTensorValue>("rock_thermal_conductivity")),
    _lambda_fluid(getParam<RealTensorValue>("fluid_thermal_conductivity")),
    _lambda_rock_qp(declareProperty<RealTensorValue>("rock_thermal_conductivity_qp")),
    _lambda_fluid_qp(declareProperty<RealTensorValue>("fluid_thermal_conductivity_qp")),
    _dlambda_rock_qp_dvar(declareProperty<std::vector<RealTensorValue>>("dlambda_rock_qp_dvar")),
    _dlambda_fluid_qp_dvar(declareProperty<std::vector<RealTensorValue>>("dlambda_fluid_qp_dvar")),
    _porosity_qp(&getMaterialProperty<Real>("PorousFlow_porosity_qp")),
    _dporosity_qp_dvar(&getMaterialProperty<std::vector<Real>>("dPorousFlow_porosity_qp_dvar")),
    _modify(getParam<bool>("modify")),
    _power(getParam<Real>("power"))
{
}

void PorousFlowThermalConductivityTensor::computeQpProperties()
{
  _lambda_rock_qp[_qp] = _lambda_rock;
  if(_modify)
    _lambda_rock_qp[_qp] = std::pow(1 - (*_porosity_qp)[_qp], _power);  

  _dlambda_rock_qp_dvar[_qp].assign(_num_var, RealTensorValue());
  if(_modify)
  {
    for(unsigned v = 0; v < _num_var; ++v)
      _dlambda_rock_qp_dvar[_qp][v] = 
         - _power * std::pow(1 - (*_porosity_qp)[_qp], _power - 1.0) *
         (*_dporosity_qp_dvar)[_qp][v];
  }

  _dlambda_fluid_qp_dvar[_qp].assign(_num_var, RealTensorValue());
  if(_modify)
  {
    for(unsigned v = 0; v < _num_var; ++v)
      _dlambda_fluid_qp_dvar[_qp][v] = 
          _power * std::pow((*_porosity_qp)[_qp], _power - 1.0) *
          (*_dporosity_qp_dvar)[_qp][v];
  }
}
