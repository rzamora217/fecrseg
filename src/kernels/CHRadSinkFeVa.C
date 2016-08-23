#include "CHRadSinkFeVa.h"
#include "MathUtils.h"
//#include "AddV.h"
using namespace MathUtils;

template<>
InputParameters validParams<CHRadSinkFeVa>()
{
   InputParameters params = validParams<Kernel>();
   
   //params.addParam<std::string>("mob_name","M","The mobility used with the kernel");
   //params.addRequiredParam<unsigned int>("crys_num", "number of grains");
   //params.addRequiredParam<unsigned int>("op_num", "number of grains, used to generate vector v");
   //params.addRequiredParam<std::string>("var_name_base","base for variable names");
   //params.addCoupledVar("v", "Array of coupled variables");
   //params.addParam<std::string>("LogC_name","LogC","The name of the Xe mobility parameter");

   params.addParam<MaterialPropertyName>("mob_name","M","The mobility used with the kernel");
   params.addRequiredCoupledVarWithAutoBuild("v", "var_name_base", "op_num", "Array of coupled variables");
   params.addParam<MaterialPropertyName>("LogC_name","LogC","The name of the Xe mobility parameter");
   params.addParam<MaterialPropertyName>("LogTol_name","LogTol","The name of the Xe mobility parameter");

   return params;
}

CHRadSinkFeVa::CHRadSinkFeVa(const InputParameters & parameters)
  :Kernel(parameters),
   _vVaC_i(getMaterialProperty<std::vector<Real> >("vVaC_i")),
   _vVam_i(getMaterialProperty<std::vector<Real> >("vVam_i")),
   _vVab_i(getMaterialProperty<std::vector<Real> >("vVab_i")),
   _mob_name(getParam<MaterialPropertyName>("mob_name")),
   _M(getMaterialProperty<Real>(_mob_name)),
   _kb(getMaterialProperty<Real>("kb")),
   _kT(getMaterialProperty<Real>("kT")),
   _LogC_name(getParam<MaterialPropertyName>("LogC_name")),
   _LogC(getMaterialProperty<Real>(_LogC_name)),
   _LogTol_name(getParam<MaterialPropertyName>("LogTol_name")),
   _LogTol(getMaterialProperty<Real>(_LogTol_name)),
   _n(coupledComponents("v"))
 
{
//  _n = coupledComponents("v");
  _vals.resize(_n);
  _grad_vals.resize(_n);
  
  for (int i=0; i<_n; ++i)
    {
      _vals[i] = &coupledValue("v", i);
      _grad_vals[i] = &coupledGradient("v", i);
    }
  
}

Real
CHRadSinkFeVa::computeQpResidual()
{
  
  Real tol = _LogTol[_qp];
  
  RealGradient outval1 = 0.0;
  
  unsigned int ind = 0;
  
  Real ytot = _u[_qp];
  RealGradient grad_ytot = _grad_u[_qp];
  
  for (int i=0; i<_n; ++i)
  {
    for (int j=i+1; j<_n; ++j)
    {
	
	  Real f_i = -8.0*(*_vals[i])[_qp]*(*_vals[i])[_qp]*(*_vals[j])[_qp]*(*_vals[j])[_qp];
	  RealGradient grad_f_i = -8.0*2.0*(*_vals[i])[_qp]*(*_grad_vals[i])[_qp]*(*_vals[j])[_qp]*(*_vals[j])[_qp];
	  grad_f_i -= 8.0*(*_vals[i])[_qp]*(*_vals[i])[_qp]*2.0*(*_vals[j])[_qp]*(*_grad_vals[j])[_qp];
      if(_vVam_i[_qp][ind]>(1e-16))
      {
	    Real inside_g_i = 100.0*ytot + exp(-_vVab_i[_qp][ind]/_vVam_i[_qp][ind]);
	    Real dgidc = 100.0*_vVam_i[_qp][ind]*poly4Log(inside_g_i,tol,1);
	    Real g_i = _vVam_i[_qp][ind]*poly4Log(inside_g_i,tol,0) + _vVab_i[_qp][ind];
	    RealGradient grad_g_i = grad_ytot*dgidc;
	    outval1 += -grad_g_i*f_i + (_vVaC_i[_qp][ind] - g_i)*grad_f_i;
      }else{
	    outval1 += _vVaC_i[_qp][ind]*grad_f_i;
      }
	  ind++;
    }
  }
  return  _grad_test[_i][_qp]*_M[_qp]*_u[_qp]*2*outval1;
  mooseError("Unsupported type passed in.");
}


Real
CHRadSinkFeVa::computeQpJacobian()
{
  
  Real tol = _LogTol[_qp];
  
  //Terms used by the Residual and Jacobian calculation

  Real ytot = _u[_qp];
  RealGradient grad_ytot= _grad_u[_qp];
  
  RealGradient outval1 = 0.0;
  RealGradient doutval1 = 0.0;
  
  unsigned int ind = 0;
  
  for (int i=0; i<_n; ++i)
  {
    for (int j=i+1; j<_n; ++j)
    {
	  Real f_i = -8.0*(*_vals[i])[_qp]*(*_vals[i])[_qp]*(*_vals[j])[_qp]*(*_vals[j])[_qp];
	  RealGradient grad_f_i = -8.0*2.0*(*_vals[i])[_qp]*(*_grad_vals[i])[_qp]*(*_vals[j])[_qp]*(*_vals[j])[_qp];
	  grad_f_i -= 8.0*(*_vals[i])[_qp]*(*_vals[i])[_qp]*2.0*(*_vals[j])[_qp]*(*_grad_vals[j])[_qp];

      if(_vVam_i[_qp][ind]>(1e-16))
      {    
	    Real inside_g_i = 100.0*ytot + exp(-_vVab_i[_qp][ind]/_vVam_i[_qp][ind]);
	    Real dgidc = 100.0*_vVam_i[_qp][ind]*poly4Log(inside_g_i,tol,1);
	    Real g_i = _vVam_i[_qp][ind]*poly4Log(inside_g_i,tol,0) + _vVab_i[_qp][ind];
	    RealGradient grad_g_i = grad_ytot*dgidc;
	    outval1 += -grad_g_i*f_i + (_vVaC_i[_qp][ind] - g_i)*grad_f_i;
	    doutval1+=(-_vVam_i[_qp][ind]*100*_phi[_j][_qp]*poly4Log((100*ytot+exp(-_vVab_i[_qp][ind]/_vVam_i[_qp][ind])),tol,1)*grad_f_i-((_vVam_i[_qp][ind]*100*_grad_phi[_j][_qp]*(100*ytot+exp(-_vVab_i[_qp][ind]/_vVam_i[_qp][ind]))-100*_phi[_j][_qp]*_vVam_i[_qp][ind]*100*grad_ytot)*poly4Log(((100*ytot+exp(-_vVab_i[_qp][ind]/_vVam_i[_qp][ind]))),tol,2)*(-1))*f_i);
      }else{
	    outval1 += _vVaC_i[_qp][ind]*grad_f_i;
      }
      ind++;
    }
  }
  return (_grad_test[_i][_qp]*_M[_qp]*_u[_qp]*2*doutval1+_grad_test[_i][_qp]*_M[_qp]*_phi[_j][_qp]*2*outval1);
  mooseError("Unsupported type passed in.");
}
