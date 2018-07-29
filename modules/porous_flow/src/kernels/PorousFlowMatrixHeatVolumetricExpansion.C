//*

#include "PorousFlowMatrixHeatVolumetricExpansion.h"

#include "MooseVariable.h"

registerMooseObject("PorousFlowApp", PorousFlowMatrixHeatVolumetricExpansion);

template <>
InputParameters
validParams<PorousFlowMatrixHeatVolumetricExpansion>()
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
  params.addClassDescription("Energy-density*rate_of_solid_volumetric_expansion");
  return params;
}

PorousFlowMatrixHeatVolumetricExpansion::PorousFlowMatrixHeatVolumetricExpansion(const InputParameters & parameters)
  : TimeKernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _var_is_porflow_var(_dictator.isPorousFlowVariable(_var.number())),
    _strain_at_nearest_qp(getParam<bool>("strain_at_nearest_qp")),
    _porosity(getMaterialProperty<Real>("PorousFlow_porosity_nodal")),
    _dporosity_dvar(getMaterialProperty<std::vector<Real>>("dPorousFlow_porosity_nodal_dvar")),
    _dporosity_dgradvar(
        getMaterialProperty<std::vector<RealGradient>>("dPorousFlow_porosity_nodal_dgradvar")),
    _nearest_qp(_strain_at_nearest_qp
                    ? &getMaterialProperty<unsigned int>("PorousFlow_nearestqp_nodal")
                    : nullptr),
    _rock_energy_nodal(getMaterialProperty<Real>("PorousFlow_matrix_internal_energy_nodal")),
    _drock_energy_nodal_dvar(
        getMaterialProperty<std::vector<Real>>("dPorousFlow_matrix_internal_energy_nodal_dvar")),
    _strain_rate_qp(getMaterialProperty<Real>("PorousFlow_volumetric_strain_rate_qp")),
    _dstrain_rate_qp_dvar(getMaterialProperty<std::vector<RealGradient>>(
        "dPorousFlow_volumetric_strain_rate_qp_dvar"))
{
}

Real
PorousFlowMatrixHeatVolumetricExpansion::computeQpResidual()
{
  Real energy = (1.0 - _porosity[_i]) * _rock_energy_nodal[_i];
  return _test[_i][_qp] * energy * _strain_rate_qp[_qp];
}

Real
PorousFlowMatrixHeatVolumetricExpansion::computeQpJacobian()
{
  return computedEnergyQpJac(_var.number()) + computedVolQpJac(_var.number());
}

Real
PorousFlowMatrixHeatVolumetricExpansion::computeQpOffDiagJacobian(unsigned int jvar)
{
  return computedEnergyQpJac(jvar) + computedVolQpJac(jvar);
}

Real
PorousFlowMatrixHeatVolumetricExpansion::computedVolQpJac(unsigned int jvar)
{
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;

  Real energy = (1.0 - _porosity[_i]) * _rock_energy_nodal[_i];

  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);
  Real dvol = _dstrain_rate_qp_dvar[_qp][pvar] * _grad_phi[_j][_qp];

  return _test[_i][_qp] * energy * dvol;
}

Real
PorousFlowMatrixHeatVolumetricExpansion::computedEnergyQpJac(unsigned int jvar)
{
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;

  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);
  const unsigned nearest_qp = (_strain_at_nearest_qp ? (*_nearest_qp)[_i] : _i);

  Real denergy = -_dporosity_dgradvar[_i][pvar] * _grad_phi[_j][_i] * _rock_energy_nodal[_i];

  if (_i != _j)
    return _test[_i][_qp] * denergy * _strain_rate_qp[_qp];

  denergy += _drock_energy_nodal_dvar[_i][pvar] * (1.0 - _porosity[_i]);
  denergy -= _rock_energy_nodal[_i] * _dporosity_dvar[_i][pvar];

  return _test[_i][_qp] * denergy * _strain_rate_qp[_qp];

}
