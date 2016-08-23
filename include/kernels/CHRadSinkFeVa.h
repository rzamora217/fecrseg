#ifndef CHRadSinkFeVa_H
#define CHRadSinkFeVa_H

#include "Kernel.h"
#include "Material.h"

//Forward Declarations
class CHRadSinkFeVa;

template<>
InputParameters validParams<CHRadSinkFeVa>();

class CHRadSinkFeVa : public Kernel
{
public:

  CHRadSinkFeVa(const InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  
  std::vector<const VariableValue *> _vals;
  std::vector<const VariableGradient *> _grad_vals;
  
  const MaterialProperty<std::vector<Real> > & _vVaC_i;
  const MaterialProperty<std::vector<Real> > & _vVam_i;
  const MaterialProperty<std::vector<Real> > & _vVab_i;

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
#endif //CHRadSinkFeVa
