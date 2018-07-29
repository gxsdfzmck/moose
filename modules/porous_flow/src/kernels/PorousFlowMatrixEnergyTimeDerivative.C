//* 

#include "PorousFlowMatrixEnergyTimeDerivative.h"

//
#include "MooseVariable.h"

registerMooseObject("PorousFlowApp", PorousFlowMatrixEnergyTimeDerivative);

template<>
InputParameters
validParams<PorousFlowMatrixEnergyTimeDerivative>()
{
  InputParameters params = validParams<TimeKernel>();
  params.addParam<bool>("strain_at_nearest_qp",
                        false,
                        "When calculating nodal porosity that depends on strain, use the strain at "
                        "the nearest quadpoint.  This adds a small extra computational burden, and "
                        "is not necessary for simulations involving only linear lagrange elements. "
                        " If you set this to true, you will also want to set the same parameter to "
                        "true for related Kernels and Materials");
  params.addRequiredParam<UserObjectName>(
      "PorousFlowDictator", "The UserObject that holds the list of Porous-Flow variable names.");
  params.addClassDescription("Derivative of heat-energy-density wrt time");
  return params;
}

PorousFlowMatrixEnergyTimeDerivative::PorousFlowMatrixEnergyTimeDerivative(const InputParameters & parameters)
  : TimeKernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _var_is_porflow_var(_dictator.isPorousFlowVariable(_var.number())),
    _strain_at_nearest_qp(getParam<bool>("strain_at_nearest_qp")),
    _porosity(getMaterialProperty<Real>("PorousFlow_porosity_nodal")),
    _porosity_old(getMaterialPropertyOld<Real>("PorousFlow_porosity_nodal")),
    _dporosity_dvar(getMaterialProperty<std::vector<Real>>("dPorousFlow_porosity_nodal_dvar")),
    _dporosity_dgradvar(
        getMaterialProperty<std::vector<RealGradient>>("dPorousFlow_porosity_nodal_dgradvar")),
    _nearest_qp(_strain_at_nearest_qp
                    ? &getMaterialProperty<unsigned int>("PorousFlow_nearestqp_nodal")
                    : nullptr),
    _rock_energy_nodal(getMaterialProperty<Real>("PorousFlow_matrix_internal_energy_nodal")),
    _rock_energy_nodal_old(getMaterialPropertyOld<Real>("PorousFlow_matrix_internal_energy_nodal")),
    _drock_energy_nodal_dvar(
        getMaterialProperty<std::vector<Real>>("dPorousFlow_matrix_internal_energy_nodal_dvar"))
{
}

Real
PorousFlowMatrixEnergyTimeDerivative::computeQpResidual()
{
  Real energy = (1.0 - _porosity[_i]) * _rock_energy_nodal[_i];
  Real energy_old = (1.0 - _porosity[_i]) * _rock_energy_nodal_old[_i];
return _test[_i][_qp] * (energy - energy_old) / _dt;
}

Real
PorousFlowMatrixEnergyTimeDerivative::computeQpJacobian()
{
  /// If the variable is not a PorousFlow variable (very unusual), the diag Jacobian terms are 0
  if (!_var_is_porflow_var)
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(_var.number()));
}

Real
PorousFlowMatrixEnergyTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  /// If the variable is not a PorousFlow variable, the OffDiag Jacobian terms are 0
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(jvar));
}

Real
PorousFlowMatrixEnergyTimeDerivative::computeQpJac(unsigned int pvar) const
{
  const unsigned nearest_qp = (_strain_at_nearest_qp ? (*_nearest_qp)[_i] : _i);

  // porosity is dependent on variables that are lumped to the nodes,
  // but it can depend on the gradient
  // of variables, which are NOT lumped to the nodes, hence:
  Real denergy = -_dporosity_dgradvar[_i][pvar] * _grad_phi[_j][_i] * _rock_energy_nodal[_i];

  if (_i != _j)
    return _test[_i][_qp] * denergy / _dt;

  /// As the fluid energy is lumped to the nodes, only non-zero terms are for _i==_j
  denergy += -_dporosity_dvar[_i][pvar] * _rock_energy_nodal[_i];
  denergy += (1.0 - _porosity[_i]) * _drock_energy_nodal_dvar[_i][pvar];

  return _test[_i][_qp] * denergy / _dt;
}
