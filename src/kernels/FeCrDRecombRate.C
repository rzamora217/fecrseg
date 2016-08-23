#include "FeCrDRecombRate.h"
#include "MathUtils.h"

template<>
InputParameters validParams<FeCrDRecombRate>()
{
  InputParameters params = validParams<Kernel>();
  params.addRequiredParam<MaterialPropertyName>("diffusivityCrV","The effective difusivity used with the kernel");
  params.addRequiredParam<MaterialPropertyName>("diffusivityFeV","The effective difusivity used with the kernel");
  params.addRequiredParam<MaterialPropertyName>("diffusivityCrI","The effective difusivity used with the kernel");
  params.addRequiredParam<MaterialPropertyName>("diffusivityFeI","The effective difusivity used with the kernel");
  params.addRequiredCoupledVar("cCr", "The coupled concentration with governing gradient");
  params.addRequiredCoupledVar("c", " Concentration");
  params.addParam<MaterialPropertyName>("LogC_name","LogC","*undocumented*");
  params.addParam<MaterialPropertyName>("LogTol_name","LogTol","*undocumented*");
  params.addRequiredCoupledVar("cVa","Concentration");
  params.addRequiredCoupledVar("cSI","Concentration");

  return params;
}

FeCrDRecombRate::FeCrDRecombRate(const InputParameters & parameters) :
    Kernel(parameters),
    _diff_nameCrV(getParam<MaterialPropertyName>("diffusivityCrV")),
    _dCrV(getMaterialProperty<Real>(_diff_nameCrV)),
    _diff_nameFeV(getParam<MaterialPropertyName>("diffusivityFeV")),
    _dFeV(getMaterialProperty<Real>(_diff_nameFeV)),
    _diff_nameCrI(getParam<MaterialPropertyName>("diffusivityCrI")),
    _dCrI(getMaterialProperty<Real>(_diff_nameCrI)),
    _diff_nameFeI(getParam<MaterialPropertyName>("diffusivityFeI")),
    _dFeI(getMaterialProperty<Real>(_diff_nameFeI)),
    _cCr_var(coupled("cCr")),
    _cCr(coupledValue("cCr")),
    _Drecomb(getMaterialProperty<Real>("Drecomb")),
    _cFe(getMaterialProperty<Real>("cFe")),
    _grad_cFe(getMaterialProperty<RealGradient>("grad_cFe")),
    _c(coupledValue("c")),
    _grad_c(coupledGradient("c")),
    _kT(getMaterialProperty<Real>("kT")),
    _Zg(getMaterialProperty<Real>("Zg")),
    _dgHVaUI(getMaterialProperty<Real>("dgHVaUI")),
    _LogC_name(getParam<MaterialPropertyName>("LogC_name")),
    _LogC(getMaterialProperty<Real>(_LogC_name)),
    _LogTol_name(getParam<MaterialPropertyName>("LogTol_name")),
    _LogTol(getMaterialProperty<Real>(_LogTol_name)),
    _kappa_cva(getMaterialProperty<Real>("kappa_cVa")),
    _kappa_csI(getMaterialProperty<Real>("kappa_cSI")),
    _kappa_ccr(getMaterialProperty<Real>("kappa_cCr")),
    _cV(coupledValue("cVa")),
    _cI(coupledValue("cSI"))
{
}

Real
FeCrDRecombRate::computeQpResidual()
{
  const Real tol = _LogTol[_qp];
  Real _S   = 0.2866; //0.250;
  Real _Vat = 0.011770599; //0.024;
  Real _Rate = 4.0*3.14159265359*_S*(_Drecomb[_qp])/_Vat;
  return 1.0*_test[_i][_qp]*(_c[_qp]*_u[_qp])*_Rate;
}

Real
FeCrDRecombRate::computeQpJacobian()
{
  const Real tol = _LogTol[_qp];
  Real _S   = 0.2866; //0.250;
  Real _Vat = 0.011770599; //0.024;
  Real _Rate = 4.0*3.14159265359*_S*(_Drecomb[_qp])/_Vat;
  return 1.0*_test[_i][_qp]*(_c[_qp]*_phi[_j][_qp])*_Rate;

}
