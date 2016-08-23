#include "FeCrIdealSink.h"
#include "MathUtils.h"
//#include "AddV.h"
using namespace MathUtils;

template<>
InputParameters validParams<FeCrIdealSink>()
{
   InputParameters params = validParams<Kernel>();
   //params.addRequiredCoupledVarWithAutoBuild("v", "var_name_base", "op_num", "Array of coupled variables");
   params.addRequiredCoupledVar("c", " Concentration");
   params.addRequiredParam<Real>("xdim", "x dimension of geometry");
   params.addRequiredParam<Real>("dose_rate", "Va and SI source rate");
   params.addRequiredParam<Real>("fcut", "cutoff for ideal sink");
   params.addRequiredParam<MaterialPropertyName>("diffusivitySink","The effective difusivity used with the kernel");
   return params;
}

FeCrIdealSink::FeCrIdealSink(const InputParameters & parameters)
  :Kernel(parameters),
   //_n(coupledComponents("v")),
   _c(coupledValue("c")),
   _xdim(getParam<Real>("xdim")),
   _dose_rate(getParam<Real>("dose_rate")),
   _fcut(getParam<Real>("fcut")),
   _diffusivitySink(getParam<MaterialPropertyName>("diffusivitySink")),
   _diffusivity(getMaterialProperty<Real>(_diffusivitySink))
{
  //_vals.resize(_n);
  //_grad_vals.resize(_n);
  //for (int i=0; i<_n; ++i)
  //{
  //  _vals[i] = &coupledValue("v", i);
  //  _grad_vals[i] = &coupledGradient("v", i);
  //}
}

Real
FeCrIdealSink::computeQpResidual()
{
  //Real ytot = _u[_qp];
  //RealGradient grad_ytot= _grad_u[_qp];
  //unsigned int ind = 0;
  
  //Real f_i = 1.0;
  //f_i -= std::abs((*_vals[0])[_qp]-(*_vals[1])[_qp]); // Assumes 2 GRAINS!!
 
  Real f_i = 0.0;
  Real fdistl = std::abs(_q_point[_qp](0)-0.0);
  Real fdistr = std::abs(_q_point[_qp](0)-_xdim);
  Real ffact  = 1.0/_fcut;
  Real _S   = 0.2866; //0.250;
  Real _Vat = 0.011770599; //0.024;
  Real _Rate = 4.0*3.14159265359*_S*_diffusivity[_qp]/_Vat;
  if(fdistl < _fcut)
  {
    //f_i = 1.0-(fdistl*ffact)*(fdistl*ffact)*(fdistl*ffact);
    f_i = 1.0;
  }
  else if(fdistr < _fcut)
  {
    //f_i = 1.0-(fdistr*ffact)*(fdistr*ffact)*(fdistr*ffact);
    f_i = 1.0;
  }
  
  if(f_i<0.0) f_i = 0.0;
  if(f_i>1.0) f_i = 1.0;
  return (_test[_i][_qp]*_u[_qp]*_Rate+_test[_i][_qp]*_dose_rate)*f_i;
  //return _test[_i][_qp]*_u[_qp]*_diffusivity[_qp]*f_i;
}

Real
FeCrIdealSink::computeQpJacobian()
{
  //Real ytot = _u[_qp];
  //RealGradient grad_ytot= _grad_u[_qp];
  //unsigned int ind = 0;
  
  //Real f_i = 1.0;
  //f_i -= std::abs((*_vals[0])[_qp]-(*_vals[1])[_qp]); // Assumes 2 GRAINS!!
 
  Real f_i = 0.0;
  Real fdistl = std::abs(_q_point[_qp](0)-0.0);
  Real fdistr = std::abs(_q_point[_qp](0)-_xdim);
  Real ffact  = 1.0/_fcut;
  Real _S   = 0.2866; //0.250;
  Real _Vat = 0.011770599; //0.024;
  Real _Rate = 4.0*3.14159265359*_S*_diffusivity[_qp]/_Vat;
  if(fdistl < _fcut)
  {
    //f_i = 1.0-(fdistl*ffact)*(fdistl*ffact)*(fdistl*ffact);
    f_i = 1.0;
  }
  else if(fdistr < _fcut)
  {
    //f_i = 1.0-(fdistr*ffact)*(fdistr*ffact)*(fdistr*ffact);
    f_i = 1.0;
  }
  
  if(f_i<0.0) f_i = 0.0;
  if(f_i>1.0) f_i = 1.0;
  return _test[_i][_qp]*_phi[_j][_qp]*_Rate*f_i;
}

