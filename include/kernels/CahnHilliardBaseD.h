/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CAHNHILLIARDBASED_H
#define CAHNHILLIARDBASED_H

#include "CHBulk.h"

/**
 * CahnHilliardBaseD implements the residual of the Cahn-Hilliard
 * equation in a general way that can be templated to a scalar or
 * tensor mobility.
 */
template<typename T>
class CahnHilliardBaseD : public CHBulk<T>
{
public:
  CahnHilliardBaseD(const InputParameters & parameters);

  static InputParameters validParams();
  virtual void initialSetup();

protected:
  typedef typename CHBulk<T>::PFFunctionType PFFunctionType;
  virtual RealGradient computeGradDFDCons(PFFunctionType type);
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  // Explicitly declare the use of the following members of the parent class
  // https://isocpp.org/wiki/faq/templates#nondependent-name-lookup-members
  using CHBulk<T>::_M;
  using CHBulk<T>::_i;
  using CHBulk<T>::_j;
  using CHBulk<T>::_qp;
  using CHBulk<T>::_var;
  using CHBulk<T>::_phi;
  using CHBulk<T>::_grad_u;
  using CHBulk<T>::_grad_phi;
  using CHBulk<T>::_grad_test;
  using CHBulk<T>::_coupled_moose_vars;

private:
  const unsigned int _nvar;
  std::vector<const MaterialProperty<Real>* > _second_derivatives;
  std::vector<const MaterialProperty<Real>* > _third_derivatives;
  std::vector<std::vector<const MaterialProperty<Real>* > > _third_cross_derivatives;
  std::vector<const VariableGradient *> _grad_vars;
};

template<typename T>
InputParameters
CahnHilliardBaseD<T>::validParams()
{
  InputParameters params = CHBulk<Real>::validParams();
  params.addClassDescription("Cahn-Hilliard Kernel that uses a DerivativeMaterial Free Energy");
  params.addRequiredParam<MaterialPropertyName>("f_name", "Base name of the free energy function F defined in a DerivativeParsedMaterial");
  params.addCoupledVar("displacement_gradients", "Vector of displacement gradient variables (see Modules/PhaseField/DisplacementGradients action)");
  return params;
}

template<typename T>
CahnHilliardBaseD<T>::CahnHilliardBaseD(const InputParameters & parameters) :
    CHBulk<T>(parameters),
    _nvar(_coupled_moose_vars.size()),
    _second_derivatives(_nvar+1),
    _third_derivatives(_nvar+1),
    _third_cross_derivatives(_nvar),
    _grad_vars(_nvar+1)
{
  // derivatives w.r.t. and gradients of the kernel variable
  _second_derivatives[0] = &this->template getMaterialPropertyDerivative<Real>("f_name", _var.name(), _var.name());
  _third_derivatives[0]  = &this->template getMaterialPropertyDerivative<Real>("f_name", _var.name(), _var.name(), _var.name());
  _grad_vars[0] = &(_grad_u);
  
  for (unsigned int i = 0; i < _nvar; ++i)
  {
    VariableName iname = _coupled_moose_vars[i]->name();
    std::cout<<"iname = "<<iname<<std::endl;
    if(iname == "cCr")
    {
      _second_derivatives[0] = &this->template getMaterialPropertyDerivative<Real>("f_name", iname, iname);
      _third_derivatives[0]  = &this->template getMaterialPropertyDerivative<Real>("f_name", iname, iname, iname);
      _grad_vars[0] = &(_coupled_moose_vars[i]->gradSln());
    }
  }
}

template<typename T>
void
CahnHilliardBaseD<T>:: initialSetup()
{
  /**
   * Check if both the non-linear as well as the auxiliary variables variables
   * are coupled. Derivatives with respect to both types of variables contribute
   * the residual.
   */
  this->template validateCoupling<Real>("f_name", "cCr"); //this->template validateCoupling<Real>("f_name", _var.name());
  //this->template validateDerivativeMaterialPropertyBase<Real>("f_name");
}

template<typename T>
RealGradient
CahnHilliardBaseD<T>::computeGradDFDCons(PFFunctionType type)
{
  RealGradient res = 0.0;
  switch (type)
  {
    case CHBulk<T>::Residual:
      res += (*_grad_vars[0])[_qp] * (*_second_derivatives[0])[_qp];
      return res;

    case CHBulk<T>::Jacobian:
      res = _grad_phi[_j][_qp] * (*_second_derivatives[0])[_qp];
      res += _phi[_j][_qp] * (*_grad_vars[0])[_qp] * (*_third_derivatives[0])[_qp];
      return res;
  }
  mooseError("Internal error");
}

template<typename T>
Real
CahnHilliardBaseD<T>::computeQpOffDiagJacobian(unsigned int jvar)
{
  return 0.0;
}

#endif // CAHNHILLIARDBASED_H
