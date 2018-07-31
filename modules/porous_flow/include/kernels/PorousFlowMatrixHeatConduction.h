//*

#ifndef MATRIXHEATCONDUCTION_H
#define MATRIXHEATCONDUCTION_H

#include "Kernel.h"
#include "PorousFlowDictator.h"

// Forward Declaractions
class PorousFlowMatrixHeatConduction;

template <>
InputParameters validParams<PorousFlowMatrixHeatConduction>();

/**
 * Kernel = grad(test) * matrix_thermal_conductivity * grad(matrix_temperatrue)
 */
class PorousFlowMatrixHeatConduction : public Kernel
{
public:
  PorousFlowMatrixHeatConduction(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  /// the Material property base name selected from: "matrix", "fluid" or "equilibrium" (default)
  std::string _base_name;

  /// holds info on the PorousFlow variables 
  const PorousFlowDictator & _dictator;
 
  /// matrix thermal conductivity at the quadpoints
  const MaterialProperty<RealTensorValue> & _lambda_rock;
 
  /// d(matrix_thermal_conductivity at qps) / d(PorousFlow variables)
  const MaterialProperty<std::vector<RealTensorValue>> & _dlambda_rock_dvar;

  /// grad(temperature)
  const MaterialProperty<RealGradient> & _grad_t;

  /// d(gradT)/d(PorousFlow variable)
  const MaterialProperty<std::vector<RealGradient>> & _dgrad_t_dvar;

  ///d(gradT)/d(grad PorousFlow variable)
  const MaterialProperty<std::vector<Real>> & _dgrad_t_dgradvar;

};

#endif
