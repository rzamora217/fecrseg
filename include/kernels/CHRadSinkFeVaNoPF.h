#ifndef CHRadSinkFeVaNoPF_H
#define CHRadSinkFeVaNoPF_H

#include "Kernel.h"
#include "Material.h"

//Forward Declarations
class CHRadSinkFeVaNoPF;

template<>
InputParameters validParams<CHRadSinkFeVaNoPF>();

class CHRadSinkFeVaNoPF : public Kernel
{
public:

  CHRadSinkFeVaNoPF(const InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  Real _xdim;
  std::string _mob_name;
  const MaterialProperty<Real> & _M;
  std::string _LogC_name;
  const MaterialProperty<Real> & _LogC;
  std::string _LogTol_name;
  const MaterialProperty<Real> & _LogTol;
  
};
#endif //CHRadSinkFeVaNoPF
