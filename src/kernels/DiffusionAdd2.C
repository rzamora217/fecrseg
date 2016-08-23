/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "DiffusionAdd2.h"

template<>
InputParameters validParams<DiffusionAdd2>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<MaterialPropertyName>("diffusivityCr","The effective difusivity used with the kernel");
  params.addRequiredParam<MaterialPropertyName>("diffusivityFe","The effective difusivity used with the kernel");
  params.addRequiredCoupledVar("cCr", "The coupled concentration with governing gradient");
  params.addRequiredCoupledVar("cVa", "The coupled concentration with governing gradient");
  params.addRequiredParam<Real>("xdim", "x dimension of geometry");
  params.addRequiredParam<Real>("fcut", "cutoff for ideal sink");
  return params;
}

DiffusionAdd2::DiffusionAdd2(const InputParameters & parameters) :
    Kernel(parameters),
    _diff_nameCr(getParam<MaterialPropertyName>("diffusivityCr")),
    _dCr(getMaterialProperty<Real>(_diff_nameCr)),
    _diff_nameFe(getParam<MaterialPropertyName>("diffusivityFe")),
    _dFe(getMaterialProperty<Real>(_diff_nameFe)),
    _cCr_var(coupled("cCr")),
    _cCr(coupledValue("cCr")),
    _cVa_var(coupled("cVa")),
    _cVa(coupledValue("cVa")),
    _xdim(getParam<Real>("xdim")),
    _fcut(getParam<Real>("fcut"))
    //_second_c(second("c")),
    //_second_u(second()),
    //_second_test(secondTest()),
    //_second_phi(secondPhi())
{
}

DiffusionAdd2::~DiffusionAdd2()
{
}

Real
DiffusionAdd2::computeQpResidual()
{
  //return (_cCr[_qp]*_dCr[_qp]+(1.0-_cCr[_qp])*_dFe[_qp])*_grad_u[_qp]*_grad_test[_i][_qp];
  return (_dCr[_qp]+_dFe[_qp])*_grad_u[_qp]*_grad_test[_i][_qp];
}

Real
DiffusionAdd2::computeQpJacobian()
{
  //return _grad_phi[_j][_qp]*(_cCr[_qp]*_dCr[_qp]+(1.0-_cCr[_qp])*_dFe[_qp])*_grad_test[_i][_qp];
  return _grad_phi[_j][_qp]*(_dCr[_qp]+_dFe[_qp])*_grad_test[_i][_qp];
}
