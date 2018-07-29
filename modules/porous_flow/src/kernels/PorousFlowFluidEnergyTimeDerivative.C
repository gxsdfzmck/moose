//* the fluid energy time derivative kernels

#include "PorousFlowFluidEnergyTimeDerivative.h"

// MOOSE includes
#include "MooseVariable.h"

registerMooseObject("PorousFlowApp", PorousFlowFluidEnergyTimeDerivative);

template <>
InputParameters
validParams<PorousFlowFluidEnergyTimeDerivative>()
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
  params.addClassDescription("Derivative of fluid heat-energy-density wrt time");
  return params;
}

PorousFlowFluidEnergyTimeDerivative::PorousFlowFluidEnergyTimeDerivative(const InputParameters & parameters)
  : TimeKernel(parameters),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _var_is_porflow_var(_dictator.isPorousFlowVariable(_var.number())),
    _num_phases(_dictator.numPhases()),
    _fluid_present(_num_phases > 0),
    _strain_at_nearest_qp(getParam<bool>("strain_at_nearest_qp")),
    _porosity(getMaterialProperty<Real>("PorousFlow_porosity_nodal")),
    _porosity_old(getMaterialPropertyOld<Real>("PorousFlow_porosity_nodal")),
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
    _fluid_density_old(
        _fluid_present
            ? &getMaterialPropertyOld<std::vector<Real>>("PorousFlow_fluid_phase_density_nodal")
            : nullptr),
    _dfluid_density_dvar(_fluid_present
                             ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                                   "dPorousFlow_fluid_phase_density_nodal_dvar")
                             : nullptr),
    _fluid_saturation_nodal(
        _fluid_present ? &getMaterialProperty<std::vector<Real>>("PorousFlow_saturation_nodal")
                       : nullptr),
    _fluid_saturation_nodal_old(
        _fluid_present ? &getMaterialPropertyOld<std::vector<Real>>("PorousFlow_saturation_nodal")
                       : nullptr),
    _dfluid_saturation_nodal_dvar(_fluid_present
                                      ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                                            "dPorousFlow_saturation_nodal_dvar")
                                      : nullptr),
    _energy_nodal(_fluid_present
                      ? &getMaterialProperty<std::vector<Real>>(
                            "PorousFlow_fluid_phase_internal_energy_nodal")
                      : nullptr),
    _energy_nodal_old(_fluid_present
                          ? &getMaterialPropertyOld<std::vector<Real>>(
                                "PorousFlow_fluid_phase_internal_energy_nodal")
                          : nullptr),
    _denergy_nodal_dvar(_fluid_present
                            ? &getMaterialProperty<std::vector<std::vector<Real>>>(
                                  "dPorousFlow_fluid_phase_internal_energy_nodal_dvar")
                            : nullptr)
{
}


Real
PorousFlowFluidEnergyTimeDerivative::computeQpResidual()
{
  Real energy = 0.0; 
  Real energy_old = 0.0; 
  for (unsigned ph = 0; ph < _num_phases; ++ph)
  {
    energy += (*_fluid_density)[_i][ph] * (*_fluid_saturation_nodal)[_i][ph] *
              (*_energy_nodal)[_i][ph] * _porosity[_i];
    energy_old += (*_fluid_density_old)[_i][ph] * (*_fluid_saturation_nodal_old)[_i][ph] *
                  (*_energy_nodal_old)[_i][ph] * _porosity_old[_i];
  }

  return _test[_i][_qp] * (energy - energy_old) / _dt;
}

Real
PorousFlowFluidEnergyTimeDerivative::computeQpJacobian()
{
  /// If the variable is not a PorousFlow variable (very unusual), the diag Jacobian terms are 0
  if (!_var_is_porflow_var)
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(_var.number()));
}

Real
PorousFlowFluidEnergyTimeDerivative::computeQpOffDiagJacobian(unsigned int jvar)
{
  /// If the variable is not a PorousFlow variable, the OffDiag Jacobian terms are 0
  if (_dictator.notPorousFlowVariable(jvar))
    return 0.0;
  return computeQpJac(_dictator.porousFlowVariableNum(jvar));
}

Real
PorousFlowFluidEnergyTimeDerivative::computeQpJac(unsigned int pvar) const
{
  const unsigned nearest_qp = (_strain_at_nearest_qp ? (*_nearest_qp)[_i] : _i);

  // porosity is dependent on variables that are lumped to the nodes,
  // but it can depend on the gradient
  // of variables, which are NOT lumped to the nodes, hence:
  Real denergy = 0.0; 
  for (unsigned ph = 0; ph < _num_phases; ++ph)
    denergy += (*_fluid_density)[_i][ph] * (*_fluid_saturation_nodal)[_i][ph] *
               (*_energy_nodal)[_i][ph] * _dporosity_dgradvar[_i][pvar] * _grad_phi[_j][nearest_qp];

  if (_i != _j)
    return _test[_i][_qp] * denergy / _dt;

  /// As the fluid energy is lumped to the nodes, only non-zero terms are for _i==_j
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
  return _test[_i][_qp] * denergy / _dt;
}
