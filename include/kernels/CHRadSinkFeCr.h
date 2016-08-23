#ifndef CHRadSinkFeCr_H
#define CHRadSinkFeCr_H

#include "Kernel.h"
#include "Material.h"

//Forward Declarations
class CHRadSinkFeCr;

template<>
InputParameters validParams<CHRadSinkFeCr>();

class CHRadSinkFeCr : public Kernel
{
public:

  CHRadSinkFeCr(const InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  
  std::vector<const VariableValue *> _vals;
  std::vector<const VariableGradient *> _grad_vals;
  
  const MaterialProperty<std::vector<Real> > & _vC_i;
  const MaterialProperty<std::vector<Real> > & _vm_i;
  const MaterialProperty<std::vector<Real> > & _vb_i;

  std::string _mob_name;

  const MaterialProperty<Real> & _M;

  const MaterialProperty<Real> & _kb;

  const MaterialProperty<Real> & _kT;

  std::string _LogC_name;
  
  const MaterialProperty<Real> & _LogC;
  
  std::string _LogTol_name;

  const MaterialProperty<Real> & _LogTol;

  unsigned int _n;
  
};
#endif //CHRadSinkFeCr
