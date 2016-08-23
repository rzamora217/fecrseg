/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef DIFFUSIONADD2_H
#define DIFFUSIONADD2_H

#include "Kernel.h"

//Forward Declarations
class DiffusionAdd2;

template<>
InputParameters validParams<DiffusionAdd2>();

class DiffusionAdd2 : public Kernel
{
public:

  DiffusionAdd2(const InputParameters & parameters);
  virtual ~DiffusionAdd2();

protected:

  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  
  std::string _diff_nameCr;
  const MaterialProperty<Real> & _dCr;
  std::string _diff_nameFe;
  const MaterialProperty<Real> & _dFe;
  unsigned int _cCr_var;
  const VariableValue & _cCr;
  unsigned int _cVa_var;
  const VariableValue & _cVa;
  Real _xdim;
  Real _fcut;

};
#endif //DIFFUSIONADD2_H
