#ifndef CHRadSinkFeSI_H
#define CHRadSinkFeSI_H

#include "Kernel.h"
#include "Material.h"

//Forward Declarations
class CHRadSinkFeSI;

template<>
InputParameters validParams<CHRadSinkFeSI>();

class CHRadSinkFeSI : public Kernel
{
public:

  CHRadSinkFeSI(const InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  
  std::vector<const VariableValue *> _vals;
  std::vector<const VariableGradient *> _grad_vals;
  
  const MaterialProperty<std::vector<Real> > & _vUIC_i;
  const MaterialProperty<std::vector<Real> > & _vUIm_i;
  const MaterialProperty<std::vector<Real> > & _vUIb_i;

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
#endif //CHRadSinkFeSI
