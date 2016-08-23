#ifndef FECRDRECOMBRATE_H
#define FECRDRECOMBRATE_H

#include "Kernel.h"

//Forward Declarations
class FeCrDRecombRate;

template<>
InputParameters validParams<FeCrDRecombRate>();

class FeCrDRecombRate : public Kernel
{
public:
  FeCrDRecombRate(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  std::string _diff_nameCrV;
  const MaterialProperty<Real> & _dCrV;
  std::string _diff_nameFeV;
  const MaterialProperty<Real> & _dFeV;
  std::string _diff_nameCrI;
  const MaterialProperty<Real> & _dCrI;
  std::string _diff_nameFeI;
  const MaterialProperty<Real> & _dFeI;
  unsigned int _cCr_var;
  const VariableValue & _cCr;
  const MaterialProperty<Real> & _Drecomb;
  const MaterialProperty<Real> & _cFe;
  const MaterialProperty<RealGradient> & _grad_cFe;
  const VariableValue & _c;
  const VariableGradient & _grad_c;
  const MaterialProperty<Real> & _kT;
  const MaterialProperty<Real> & _Zg;
  const MaterialProperty<Real> & _dgHVaUI;
  std::string _LogC_name;
  const MaterialProperty<Real> & _LogC;
  std::string _LogTol_name;
  const MaterialProperty<Real> & _LogTol;
  const MaterialProperty<Real> & _kappa_cva;
  const MaterialProperty<Real> & _kappa_csI;
  const MaterialProperty<Real> & _kappa_ccr;
  const VariableValue & _cV;
  const VariableValue & _cI;
};

#endif //FECRDRECOMBRATE_H
