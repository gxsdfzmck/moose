//*

#include "PorousFlowFluidHeatVolumetricExpansion.h"

#include "MooseVariable.h"

registerMooseObject("PorousFlowApp", PorousFlowFluidHeatVolumetricExpansion);

template <>
InputParameters
validParams<PorousFlowFluidHeatVolumetricExpansion>()
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
  params.addClassDescription("Fluid-energy-density*rate_of_solid_volumetric_expansion");
  return params;
}

PorousFlowFluidHeatVolumetricExpansion::PorousFlowFluidHeatVolumetricExpansion(const InputParameters & parameters)
  : TimeKernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _var_is_porflow_var(_dictator.isPorousFlowVariable(_var.number())),
    _num_phases(_dictator.numPhases()),
    _fluid_present(_num_phases > 0),
    _strain_at_nearest_qp(getParam<bool>("strain_at_nearest_qp")),
    _porosity(getMaterialProperty<Real>("PorousFlow_porosity_nodal")),
    _dporosity_dvar(getMaterialProperty<std::vector<Real>>("dPorousFlow_porosity_nodal_dvar")),
    _dporosity_dgradvar(
        getMaterialProperty<std::vector<RealGradient>>("dPorousFlow_porosity_nodal_dgradvar")),
    _nearest_qp(_strain_at_nearest_qp
                    ? &getMaterialProperty<unsigned int>("PorousFlow_nearestqp_nodal")
                    : nullptr),
    _fluid_density(
        _fluid_present
            ? &getMaterialProperty<std::vector<Real>>("PorousFlow_fluid_phase_density_nodal")
            : nullptr),
    _dfluid_density_dvar(_fluid_present
                             ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                                   "dPorousFlow_fluid_phase_density_nodal_dvar")
                             : nullptr),
    _fluid_saturation_nodal(
        _fluid_present ? &getMaterialProperty<std::vector<Real>>("PorousFlow_saturation_nodal")
                       : nullptr),
    _dfluid_saturation_nodal_dvar(_fluid_present
                                      ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                                            "dPorousFlow_saturation_nodal_dvar")
                                      : nullptr),
    _energy_nodal(_fluid_present
                      ? &getMaterialProperty<std::vector<Real>>(
                            "PorousFlow_fluid_phase_internal_energy_nodal")
                      : nullptr),
    _denergy_nodal_dvar(_fluid_present
                            ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                                  "dPorousFlow_fluid_phase_internal_energy_nodal_dvar")
                            : nullptr),
    _strain_rate_qp(getMaterialProperty<Real>("PorousFlow_volumetric_strain_rate_qp")),
    _dstrain_rate_qp_dvar(getMaterialProperty<std::vector<RealGradient>>(
        "dPorousFlow_volumetric_strain_rate_qp_dvar"))
{
}

Real
PorousFlowFluidHeatVolumetricExpansion::computeQpResidual()
{
  Real energy = 0.0;
  for (unsigned ph = 0; ph < _num_phases; ++ph)
    energy += (*_fluid_density)[_i][ph] * (*_fluid_saturation_nodal)[_i][ph] *
              (*_energy_nodal)[_i][ph] * _porosity[_i];

  return _test[_i][_qp] * energy * _strain_rate_qp[_qp];
}

Real
PorousFlowFluidHeatVolumetricExpansion::computeQpJacobian()
{
  return computedEnergyQpJac(_var.number()) + computedVolQpJac(_var.number());

}

Real
PorousFlowFluidHeatVolumetricExpansion::computeQpOffDiagJacobian(unsigned int jvar)
{
  return computedEnergyQpJac(jvar) + computedVolQpJac(jvar);

}

Real
PorousFlowFluidHeatVolumetricExpansion::computedVolQpJac(unsigned int jvar)
{
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  
  Real energy = 0.0;
    for (unsigned ph = 0; ph < _num_phases; ++ph)
    energy += (*_fluid_density)[_i][ph] * (*_fluid_saturation_nodal)[_i][ph] *
              (*_energy_nodal)[_i][ph] * _porosity[_i];

  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);
  Real dvol = _dstrain_rate_qp_dvar[_qp][pvar] * _grad_phi[_j][_qp];
  
  return _test[_i][_qp] * energy * dvol;
}

Real
PorousFlowFluidHeatVolumetricExpansion::computedEnergyQpJac(unsigned int jvar)
{
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;

  const unsigned int pvar = _dictator.porousFlowVariableNum(jvar);
  const unsigned nearest_qp = (_strain_at_nearest_qp ? (*_nearest_qp)[_i] : _i);

  Real denergy = 0.0;
  for (unsigned ph = 0; ph < _num_phases; ++ph)
  denergy += (*_fluid_density)[_i][ph] * (*_fluid_saturation_nodal)[_i][ph] *
             (*_energy_nodal)[_i][ph] * _dporosity_dgradvar[_i][pvar] * _grad_phi[_j][nearest_qp];

  if (_i != _j)
    return _test[_i][_qp] * denergy * _strain_rate_qp[_qp];

  for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
    denergy += (*_dfluid_density_dvar)[_i][ph][pvar] * (*_fluid_saturation_nodal)[_i][ph] *
               (*_energy_nodal)[_i][ph] * _porosity[_i];
    denergy += (*_fluid_density)[_i][ph] * (*_dfluid_saturation_nodal_dvar)[_i][ph][pvar] *
               (*_energy_nodal)[_i][ph] * _porosity[_i];
    denergy += (*_fluid_density)[_i][ph] * (*_fluid_saturation_nodal)[_i][ph] *
               (*_denergy_nodal_dvar)[_i][ph][pvar] * _porosity[_i];
    denergy += (*_fluid_density)[_i][ph] * (*_fluid_saturation_nodal)[_i][ph] *
               (*_energy_nodal)[_i][ph] * _dporosity_dvar[_i][pvar];
  }

  return _test[_i][_qp] * denergy * _strain_rate_qp[_qp];
}
