#ifndef FeCrVaSIbulk_H
#define FeCrVaSIbulk_H

#include "Material.h"

//Forward Declarations
class FeCrVaSIbulk;

template<>
InputParameters validParams<FeCrVaSIbulk>();

class FeCrVaSIbulk : public Material
{
public:
  FeCrVaSIbulk(const InputParameters & parameters);
  
protected:
  virtual void computeProperties();

private:
  std::vector<const VariableValue *> _vals;
  std::vector<const VariableGradient *> _grad_vals;
  std::vector<const VariableValue *> _valsVa;
  std::vector<const VariableGradient *> _grad_valsVa;

  const VariableValue & _cCr;
  const VariableValue & _cVa;
  const VariableValue & _cSI;
  const VariableGradient & _grad_cCr;
  const VariableGradient & _grad_cVa;
  const VariableGradient & _grad_cSI;
  Real _temp;
  Real _kappa_cfe;
  Real _kappa_cva;
  Real _kappa_csI;
  Real _kappa_ccr;
  Real _ICrva;
  Real _ICrsI;
  Real _IVava;
  Real _Laav1;
  Real _Lbbv1;
  Real _Labv1;
  Real _Laai1;
  Real _Lbbi1;
  Real _Labi1;
  Real _Z;
  Real _dHVaUI;
  MaterialProperty<Real> & _DrecombVa;
  MaterialProperty<Real> & _DrecombSI;
  MaterialProperty<Real> & _Drecomb;
  MaterialProperty<Real> & _dCrV;
  MaterialProperty<Real> & _dCrI;
  MaterialProperty<Real> & _dFeV;
  MaterialProperty<Real> & _dFeI;
  MaterialProperty<Real> & _dCrV_c;
  MaterialProperty<Real> & _dCrI_c;
  MaterialProperty<Real> & _dFeV_c;
  MaterialProperty<Real> & _dFeI_c;
  MaterialProperty<Real> & _kT;
  MaterialProperty<Real> & _Temp;
  MaterialProperty<Real> & _Zg;
  MaterialProperty<Real> & _dgHVaUI;
  MaterialProperty<Real> & _M_Cr;
  MaterialProperty<Real> & _M_Va;
  MaterialProperty<Real> & _M_SI;
  MaterialProperty<Real> & _M_Va_sink;
  MaterialProperty<Real> & _M_SI_sink;
  MaterialProperty<Real> & _M_Va_int;
  MaterialProperty<Real> & _M_SI_int;
  MaterialProperty<Real> & _DM_Cr;
  MaterialProperty<Real> & _DM_Va;
  MaterialProperty<Real> & _DM_SI;
  MaterialProperty<RealGradient> & _Dgrad_Mnp_Cr;
  MaterialProperty<RealGradient> & _Dgrad_Mnp_Va;
  MaterialProperty<RealGradient> & _Dgrad_Mnp_SI;
  MaterialProperty<Real> & _Dgrad_Mngp_Cr;
  MaterialProperty<Real> & _Dgrad_Mngp_Va;
  MaterialProperty<Real> & _Dgrad_Mngp_SI;
  MaterialProperty<RealGradient> & _grad_M_Cr;
  MaterialProperty<RealGradient> & _grad_M_Va;
  MaterialProperty<RealGradient> & _grad_M_SI;
  MaterialProperty<Real> & _kappa_cFe;
  MaterialProperty<Real> & _kappa_cCr;
  MaterialProperty<Real> & _kappa_cVa;
  MaterialProperty<Real> & _kappa_cSI;
  MaterialProperty<Real> & _cFe;
  MaterialProperty<RealGradient> & _grad_cFe;
  MaterialProperty<Real> & _ICrVa;
  MaterialProperty<Real> & _ICrSI;
  MaterialProperty<Real> & _IVaVa;
  MaterialProperty<Real> & _LogC;
  MaterialProperty<Real> & _LogC2;
  MaterialProperty<Real> & _LogTol;
  MaterialProperty<Real> & _LNcVa;
  MaterialProperty<Real> & _LNcSI;
  MaterialProperty<RealGradient> & _grad_LNcVa;
  MaterialProperty<RealGradient> & _grad_LNcSI;
};

#endif //FeCrVaSIbulk_H
