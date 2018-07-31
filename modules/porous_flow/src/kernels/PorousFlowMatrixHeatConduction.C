//*

#include "PorousFlowMatrixHeatConduction.h"

#include "MooseVariable.h"

registerMooseObject("PorousFlowApp", PorousFlowMatrixHeatConduction);

template<>
InputParameters
validParams<PorousFlowMatrixHeatConduction>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<UserObjectName>("PorousFlowDictator", "the UserObject that holds the list of PorousFlow variable names");
  params.addParam<std::string>("base_name", "Material property base name");
  params.addClassDescription("Matrix heat conduction in the porous Flow module");
  return params;
}

PorousFlowMatrixHeatConduction::PorousFlowMatrixHeatConduction(const InputParameters & parameters)
  : Kernel(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _lambda_rock(getMaterialProperty<RealTensorValue>("rock_thermal_conductivity_qp")),
    _dlambda_rock_dvar(getMaterialProperty<std::vector<RealTensorValue>>("dlambda_rock_qp_dvar")),
    _grad_t(getMaterialProperty<RealGradient>(_base_name + "PorousFlow_grad_temperature_qp")),
    _dgrad_t_dvar(getMaterialProperty<std::vector<RealGradient>>(_base_name + "dPorousFlow_grad_temperature_qp_dvar")),
    _dgrad_t_dgradvar(getMaterialProperty<std::vector<Real>>(_base_name + "dPorousFlow_grad_temperature_qp_dgradvar"))
{
}

Real
PorousFlowMatrixHeatConduction::computeQpResidual()
{
  return _grad_test[_i][_qp] * (_lambda_rock[_qp] * _grad_t[_qp]);
}

Real
PorousFlowMatrixHeatConduction::computeQpJacobian()
{
  return computeQpOffDiagJacobian(_var.number());
}

Real
PorousFlowMatrixHeatConduction::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;

  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);

  return  _grad_test[_i][_qp] *
         ((_dlambda_rock_dvar[_qp][pvar] * _grad_t[_qp] + _lambda_rock[_qp] * _dgrad_t_dvar[_qp][pvar]) *
              _phi[_j][_qp] +
          _lambda_rock[_qp] * _dgrad_t_dgradvar[_qp][pvar] * _grad_phi[_j][_qp]);
}
