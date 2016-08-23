#include "CHRadSinkFeSINoPF.h"
#include "MathUtils.h"
//#include "AddV.h"
using namespace MathUtils;

template<>
InputParameters validParams<CHRadSinkFeSINoPF>()
{
   InputParameters params = validParams<Kernel>();
   params.addRequiredParam<Real>("xdim", "x dimension of geometry");
   params.addParam<MaterialPropertyName>("mob_name","M","The mobility used with the kernel");
   params.addParam<MaterialPropertyName>("LogC_name","LogC","The name of the Xe mobility parameter");
   params.addParam<MaterialPropertyName>("LogTol_name","LogTol","The name of the Xe mobility parameter");

   return params;
}

CHRadSinkFeSINoPF::CHRadSinkFeSINoPF(const InputParameters & parameters) :
   Kernel(parameters),
   _xdim(getParam<Real>("xdim")),
   _mob_name(getParam<MaterialPropertyName>("mob_name")),
   _M(getMaterialProperty<Real>(_mob_name)),
   _LogC_name(getParam<MaterialPropertyName>("LogC_name")),
   _LogC(getMaterialProperty<Real>(_LogC_name)),
   _LogTol_name(getParam<MaterialPropertyName>("LogTol_name")),
   _LogTol(getMaterialProperty<Real>(_LogTol_name))
{
}

Real
CHRadSinkFeSINoPF::computeQpResidual()
{
  Real tol = _LogTol[_qp];
  RealGradient outval1 = 0.0;
  Real ytot = _u[_qp];
  RealGradient grad_ytot = _grad_u[_qp];
  
  // SI - Sigma3 GB:
  Real _C_param =  3.013 ;
  Real _k_param =  4.000 ;
  Real _b_param = -8.200 ;
  Real _m_param =  2.100 ;
  
  // Get _dist:
  Real _dist;
  Real fdistl = std::abs(_q_point[_qp](0)-0.0);
  Real fdistr = std::abs(_q_point[_qp](0)-_xdim);
  if(fdistl < fdistr){
    _dist = fdistl;
  }else{
    _dist = fdistr;
  }
  Real f_i = 1.0/(1.0+std::exp(-_dist/(_k_param*_k_param)))-1.0;
  RealGradient grad_f_i = 2.0*std::pow(1.0+std::exp(-_dist*_dist/(_k_param*_k_param)),-2.0)*std::exp(-_dist*_dist/(_k_param*_k_param))*_dist/_k_param;
  Real inside_g_i = 100.0*ytot + exp(-_b_param/_m_param);
  Real dgidc = 100.0*_m_param*poly4Log(inside_g_i,tol,1);
  Real g_i = _m_param*poly4Log(inside_g_i,tol,0) + _b_param;
  RealGradient grad_g_i = grad_ytot*dgidc;
  outval1 += -grad_g_i*f_i + (_C_param - g_i)*grad_f_i;
  return  _grad_test[_i][_qp]*_M[_qp]*_u[_qp]*2.0*outval1;
  mooseError("Unsupported type passed in.");
}


Real
CHRadSinkFeSINoPF::computeQpJacobian()
{
  Real tol = _LogTol[_qp];
  RealGradient outval1 = 0.0;
  Real ytot = _u[_qp];
  RealGradient grad_ytot = _grad_u[_qp];
  
  // SI - Sigma3 GB:
  Real _C_param =  3.013 ;
  Real _k_param =  4.000 ;
  Real _b_param = -8.200 ;
  Real _m_param =  2.100 ;
  
  // Get _dist:
  Real _dist;
  Real fdistl = std::abs(_q_point[_qp](0)-0.0);
  Real fdistr = std::abs(_q_point[_qp](0)-_xdim);
  if(fdistl < fdistr){
    _dist = fdistl;
  }else{
    _dist = fdistr;
  }
  Real f_i = 1.0/(1.0+std::exp(-_dist/(_k_param*_k_param)))-1.0;
  RealGradient grad_f_i = 2.0*std::pow(1.0+std::exp(-_dist*_dist/(_k_param*_k_param)),-2.0)*std::exp(-_dist*_dist/(_k_param*_k_param))*_dist/_k_param;
  Real inside_g_i = 100.0*ytot + exp(-_b_param/_m_param);
  Real dgidc = 100.0*_m_param*poly4Log(inside_g_i,tol,1);
  Real g_i = _m_param*poly4Log(inside_g_i,tol,0) + _b_param;
  RealGradient grad_g_i = grad_ytot*dgidc;
  outval1 += -grad_g_i*f_i + (_C_param - g_i)*grad_f_i;
  return  _grad_test[_i][_qp]*_M[_qp]*2.0*outval1;
  mooseError("Unsupported type passed in.");
}
