#ifndef POROUSFLOWTHERMALCONDUCTIVITYTENSOR_H
#define POROUSFLOWTHERMALCONDUCTIVITYTENSOR_H

#include "PorousFlowMaterialVectorBase.h"

class PorousFlowThermalConductivityTensor;

template <>
InputParameters validParams<PorousFlowThermalConductivityTensor>();

/**
 * provide the thermal conductivity tensor for the matrix and fluid 
 */
class PorousFlowThermalConductivityTensor : public PorousFlowMaterialVectorBase
{
public:
  PorousFlowThermalConductivityTensor(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  /// thermal conductivity of rock
  const RealTensorValue _lambda_rock;
 
  /// thermal conductivity of fluid
  const RealTensorValue _lambda_fluid;
  
  /// thermal conductivity of rock at the qps 
  MaterialProperty<RealTensorValue> & _lambda_rock_qp;

  /// thermal conductivity of fluid at the qps 
  MaterialProperty<RealTensorValue> & _lambda_fluid_qp;

  /// d(thermal conductivity of rock at the qps)/d(PorousFlow variable)
  MaterialProperty<std::vector<RealTensorValue>> & _dlambda_rock_qp_dvar;

  /// d(thermal conductivity of fluid at the qps)/d(PorousFlow variable)
  MaterialProperty<std::vector<RealTensorValue>> & _dlambda_fluid_qp_dvar;
 
  /// porosity of the porous media at qps
  const MaterialProperty<Real> * const _porosity_qp;

  /// d(porosity)/d(PorousFlow variable)
  const MaterialProperty<std::vector<Real>> * const _dporosity_qp_dvar;

  /// whether modify the thermal conductivity with Bruggeman correction
  const bool _modify;

  /// the power component of the correction
  const Real  _power;
};

#endif
