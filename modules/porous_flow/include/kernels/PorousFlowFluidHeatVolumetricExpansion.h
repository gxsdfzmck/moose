//*

#ifndef FLUIDHEATVOLUMETRICEXPANSION_H
#define FLUIDHEATVOLUMETRICEXPANSION_H

#include "TimeDerivative.h"
#include "PorousFlowDictator.h"

// Forward Declaration
class PorousFlowFluidHeatVolumetricExpansion;

template<>
InputParameters validParams<PorousFlowFluidHeatVolumetricExpansion>();

/**
 * Kernel = fluid_energy_density * d(volumetric_strain)/dt
 */
class PorousFlowFluidHeatVolumetricExpansion : public TimeKernel
{
public:
  PorousFlowFluidHeatVolumetricExpansion(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// holds info on the Porous Flow variables
  const PorousFlowDictator & _dictator;

  /// whether the Variable for this Kernel is a porous-flow variable according to the Dictator
  const bool _var_is_porflow_var;

  /// number of fluid phases
  const unsigned int _num_phases;

  /// whether fluid is present
  const bool _fluid_present;

  /// whether the porosity uses the volumetric strain at the closest quadpoint
  const bool _strain_at_nearest_qp;

  /// porosity
  const MaterialProperty<Real> & _porosity;

  /// d(porosity)/d(porous-flow variable)
  const MaterialProperty<std::vector<Real>> & _dporosity_dvar;

  /// d(porosity)/d(grad porous-flow variable)
  const MaterialProperty<std::vector<RealGradient>> & _dporosity_dgradvar;

  /// the nearest qp to the node
  const MaterialProperty<unsigned int> * const _nearest_qp;

  /// nodal fluid density
  const MaterialProperty<std::vector<Real>> * const _fluid_density;

  /// d(nodal fluid density)/d(porous-flow variable)
  const MaterialProperty<std::vector<std::vector<Real>>> * const _dfluid_density_dvar;

  /// nodal fluid saturation
  const MaterialProperty<std::vector<Real>> * const _fluid_saturation_nodal;

  /// d(nodal fluid saturation)/d(porous-flow variable)
  const MaterialProperty<std::vector<std::vector<Real>>> * const _dfluid_saturation_nodal_dvar;

  /// internal energy of the phases, evaluated at the nodes
  const MaterialProperty<std::vector<Real>> * const _energy_nodal;

  /// d(internal energy)/d(PorousFlow variable)
  const MaterialProperty<std::vector<std::vector<Real>>> * const _denergy_nodal_dvar;

  /// strain rate
  const MaterialProperty<Real> & _strain_rate_qp;

  /// d(strain rate)/d(porous-flow variable)
  const MaterialProperty<std::vector<RealGradient>> & _dstrain_rate_qp_dvar;

  /**
   * Derivative of energy part of the residual with respect to the Variable
   * with variable number jvar.
   * This is used by both computeQpJacobian and computeQpOffDiagJacobian
   * @param jvar take the derivative of the energy part of the residual wrt this variable number
   */
  Real computedEnergyQpJac(unsigned int jvar);

  /**
   * Derivative of volumetric-strain part of the residual with respect to the Variable
   * with variable number jvar.
   * This is used by both computeQpJacobian and computeQpOffDiagJacobian
   * @param jvar take the derivative of the volumetric-strain part of the residual wrt this variable
   * number
   */
  Real computedVolQpJac(unsigned int jvar);
};

#endif
