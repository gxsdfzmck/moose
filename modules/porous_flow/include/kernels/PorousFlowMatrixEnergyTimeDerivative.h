//*

#ifndef POROUSFLOWMATRIXENERGYTIMEDERIVATIVE
#define POROUSFLOWMATRIXENERGYTIMEDERIVATIVE

#include "TimeDerivative.h"
#include "PorousFlowDictator.h"

// Forward Delcaration

class PorousFlowMatrixEnergyTimeDerivative;

template <>
InputParameters validParams<PorousFlowMatrixEnergyTimeDerivative>();

/**
 * Kernel = (rock_heat_energy - rock_heat_energy_old)/dt
 */
class PorousFlowMatrixEnergyTimeDerivative : public TimeKernel
{
public:
  PorousFlowMatrixEnergyTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// holds info on the PorousFlow variables
  const PorousFlowDictator & _dictator;

  /// whether the Variable for this Kernel is a porous-flow variable according to the Dictator
  const bool _var_is_porflow_var;

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

  /// nodal rock energy density
  const MaterialProperty<Real> & _rock_energy_nodal;

  /// old value of nodal rock energy density
  const MaterialProperty<Real> & _rock_energy_nodal_old;

 /// d(nodal rock energy density)/d(PorousFlow variable)
  const MaterialProperty<std::vector<Real>> & _drock_energy_nodal_dvar;

 /**
   * Derivative of residual with respect to PorousFlow variable number pvar
   * This is used by both computeQpJacobian and computeQpOffDiagJacobian
   * @param pvar take the derivative of the residual wrt this PorousFlow variable
   */
  Real computeQpJac(unsigned int pvar) const;
};

#endif
