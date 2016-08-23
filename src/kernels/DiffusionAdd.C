/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "DiffusionAdd.h"

template<>
InputParameters validParams<DiffusionAdd>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<MaterialPropertyName>("diffusivity","The effective difusivity used with the kernel");
  params.addRequiredCoupledVar("c", "The coupled concentration with governing gradient");
  params.addRequiredParam<Real>("dirscale", "sign of kernel");
  return params;
}

DiffusionAdd::DiffusionAdd(const InputParameters & parameters) :
    Kernel(parameters),
    _diff_name(getParam<MaterialPropertyName>("diffusivity")),
    _diffusivity(getMaterialProperty<Real>(_diff_name)),
    _c_var(coupled("c")),
    _c(coupledValue("c")),
    _grad_c(coupledGradient("c")),
    _dirscale(getParam<Real>("dirscale"))
    //_second_c(second("c")),
    //_second_u(second()),
    //_second_test(secondTest()),
    //_second_phi(secondPhi())
{
}

DiffusionAdd::~DiffusionAdd()
{
}

Real
DiffusionAdd::computeQpResidual()
{
  // -(grad_test_i, u * d * grad_c)
  //return _dirscale*_grad_test[_i][_qp] * _u[_qp]*_diffusivity[_qp]*_grad_c[_qp];
  return _dirscale*_grad_test[_i][_qp]*_diffusivity[_qp]*_grad_c[_qp];
}

Real
DiffusionAdd::computeQpJacobian()
{
  Real value = 0.0;
  //value += _dirscale*_grad_test[_i][_qp] * _phi[_j][_qp]*_diffusivity[_qp]*_grad_c[_qp];
  value += _dirscale*_grad_test[_i][_qp]*_diffusivity[_qp]*_grad_c[_qp];
  return value;
}

