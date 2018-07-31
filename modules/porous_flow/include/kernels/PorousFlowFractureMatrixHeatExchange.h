//* This file is part of the geothermal Applicatin

#ifndef FractureMatrixHeatExchange_H
#define FractureMatrixHeatExchange_H

#include "Kernel.h"

class PorousFlowFractureMatrixHeatExchange;

template <>
InputParameters validParams<PorousFlowFractureMatrixHeatExchange>();

class PorousFlowFractureMatrixHeatExchange : public Kernel
{
public:
  PorousFlowFractureMatrixHeatExchange(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;

  const enum class HeatChangeType {
    FluidToMatrix,
    MatrixToFluid,
  } _heat_change_type;
  
private: 
  Real _h;
  const VariableValue & _couple_temperature;
};

#endif // FractureMatrixHeatExchange_H
