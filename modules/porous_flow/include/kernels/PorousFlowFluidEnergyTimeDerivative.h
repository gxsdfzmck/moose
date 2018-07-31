#ifndef POROUSFLOWFLUIDENERGYTIMEDERIVATIVE_H
#define POROUSFLOWFLUIDENERGYTIMEDERIVATIVE_H

#include "TimeDerivative.h"
#include "PorousFlowDictator.h"

// Forward Declarations
class PorousFlowFluidEnergyTimeDerivative;

template <>
InputParameters validParams<PorousFlowFluidEnergyTimeDerivative>();

/**
 * Kernel = (fluid_heat_engergy - fluid_heat_energy)/dt
 */
class PorousFlowFluidEnergyTimeDerivative : public TimeKernel
{
public:
  PorousFlowFluidEnergyTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jar) override;

  /// holds info on the PorousFlow variables
  const PorousFlowDictator & _dictator;

  /// whether the Variable for this Kernel is a porous-flow variable according to the Dictator
  const bool _var_is_porflow_var;

  /// number of fluid phases
  const unsigned int _num_phases;

  /// whether _num_phases > 0
  const bool _fluid_present;

  /// whether the porosity uses the volumetric strain at the closest quadpoint
  const bool _strain_at_nearest_qp;

  /// porosity at the nodes, but it can depend on grad(variables) which are actually evaluated at the qps
  const MaterialProperty<Real> & _porosity;

  /// old value of porosity
  const MaterialProperty<Real> & _porosity_old;

  /// d(porosity)/d(porous-flow variable) - these derivatives will be wrt variables at the nodes
  const MaterialProperty<std::vector<Real>> & _dporosity_dvar;

  /// d(porosity)/d(grad porous-flow variable) - remember these derivatives will be wrt grad(vars) at qps
  const MaterialProperty<std::vector<RealGradient>> & _dporosity_dgradvar;

  /// the nearest qp to the node
  const MaterialProperty<unsigned int> * const _nearest_qp;

  /// nodal fluid density
  const MaterialProperty<std::vector<Real>> * const _fluid_density;

  /// old value of nodal fluid density
  const MaterialProperty<std::vector<Real>> * const _fluid_density_old;

  /// d(nodal fluid density)/d(porous-flow variable)
  const MaterialProperty<std::vector<std::vector<Real>>> * const _dfluid_density_dvar;

  /// nodal fluid saturation
  const MaterialProperty<std::vector<Real>> * const _fluid_saturation_nodal;

  /// old value of fluid saturation
  const MaterialProperty<std::vector<Real>> * const _fluid_saturation_nodal_old;

  /// d(nodal fluid saturation)/d(porous-flow variable)
  const MaterialProperty<std::vector<std::vector<Real>>> * const _dfluid_saturation_nodal_dvar;

  /// internal energy of the phases, evaluated at the nodes
  const MaterialProperty<std::vector<Real>> * const _energy_nodal;

  /// old value of internal energy of the phases, evaluated at the nodes
  const MaterialProperty<std::vector<Real>> * const _energy_nodal_old;

  /// d(internal energy)/d(PorousFlow variable)
  const MaterialProperty<std::vector<std::vector<Real>>> * const _denergy_nodal_dvar;

  /**
   * Derivative of residual with respect to PorousFlow variable number pvar
   * This is used by both computeQpJacobian and computeQpOffDiagJacobian
   * @param pvar take the derivative of the residual wrt this PorousFlow variable
   */
  Real computeQpJac(unsigned int pvar) const;
};

#endif // POROUSFLOWFLUIDENERGYTIMEDERIVATIVE_H
