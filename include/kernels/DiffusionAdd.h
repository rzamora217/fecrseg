/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef DIFFUSIONADD_H
#define DIFFUSIONADD_H

#include "Kernel.h"

//Forward Declarations
class DiffusionAdd;

template<>
InputParameters validParams<DiffusionAdd>();

class DiffusionAdd : public Kernel
{
public:

  DiffusionAdd(const InputParameters & parameters);
  virtual ~DiffusionAdd();

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  
  std::string _diff_name;
  const MaterialProperty<Real> & _diffusivity;
  unsigned int _c_var;
  const VariableValue & _c;
  const VariableGradient & _grad_c;
  Real _dirscale;
  //VariableSecond & _second_c;
  //VariableSecond & _second_u;
  //VariableTestSecond & _second_test;
  //VariablePhiSecond & _second_phi;

};
#endif //DIFFUSIONADD_H
