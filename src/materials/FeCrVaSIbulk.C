#include "FeCrVaSIbulk.h"
#include "MathUtils.h"
// libMesh includes
#include "libmesh/quadrature.h"
using namespace MathUtils;

template<>
InputParameters validParams<FeCrVaSIbulk>()
{
  InputParameters params = validParams<Material>();
  params.addRequiredParam<Real>("temp", "Constant temperature in Kelvin");
  params.addParam<Real>("kappa_cfe",1.0, "The kappa parameter for the iron concentration");
  params.addParam<Real>("kappa_cva",1.0, "The kappa parameter for the vacancy concentration");
  params.addParam<Real>("kappa_ccr",1.0, "The kappa parameter for the Chromium concentration");
  params.addParam<Real>("kappa_csI",1.0, "The kappa parameter for the interstitial concentration");
  params.addParam<Real>("ICrva",0.0,"Regular solution parameter T in kB");
  params.addParam<Real>("ICrsI",0.0,"Regular solution parameter T in kB");
  params.addParam<Real>("IVava",0.0,"Regular solution parameter T in kB");
  params.addRequiredParam<Real>("Laav1","Onsager Coefficient");
  params.addRequiredParam<Real>("Lbbv1","Onsager Coefficient");
  params.addRequiredParam<Real>("Labv1","Onsager Coefficient");
  params.addRequiredParam<Real>("Laai1","Onsager Coefficient");
  params.addRequiredParam<Real>("Lbbi1","Onsager Coefficient");
  params.addRequiredParam<Real>("Labi1","Onsager Coefficient");
  params.addRequiredParam<Real>("Z", "Recombination parameter");
  params.addRequiredParam<Real>("dHVaUI", "Recombination energy");
  params.addRequiredCoupledVar("cCr","Concentration");
  params.addRequiredCoupledVar("cVa","Concentration");
  params.addRequiredCoupledVar("cSI","Concentration");
  return params;
}

FeCrVaSIbulk::FeCrVaSIbulk(const InputParameters & parameters)
  :Material(parameters),
   _cCr(coupledValue("cCr")),
   _cVa(coupledValue("cVa")),
   _cSI(coupledValue("cSI")),
   _grad_cCr(coupledGradient("cCr")),
   _grad_cVa(coupledGradient("cVa")),
   _grad_cSI(coupledGradient("cSI")),
   _temp(getParam<Real>("temp")),
   _kappa_cfe(getParam<Real>("kappa_cfe")),
   _kappa_cva(getParam<Real>("kappa_cva")),
   _kappa_csI(getParam<Real>("kappa_csI")),
   _kappa_ccr(getParam<Real>("kappa_ccr")),
   _ICrva(getParam<Real>("ICrva")),
   _ICrsI(getParam<Real>("ICrsI")),
   _IVava(getParam<Real>("IVava")),
   _Laav1(getParam<Real>("Laav1")),
   _Lbbv1(getParam<Real>("Lbbv1")),
   _Labv1(getParam<Real>("Labv1")),
   _Laai1(getParam<Real>("Laai1")),
   _Lbbi1(getParam<Real>("Lbbi1")),
   _Labi1(getParam<Real>("Labi1")),
   _Z(getParam<Real>("Z")),
   _dHVaUI(getParam<Real>("dHVaUI")),
   _DrecombVa(declareProperty<Real>("DrecombVa")),
   _DrecombSI(declareProperty<Real>("DrecombSI")),
   _Drecomb(declareProperty<Real>("Drecomb")),
   _dCrV(declareProperty<Real>("dCrV")),
   _dCrI(declareProperty<Real>("dCrI")),
   _dFeV(declareProperty<Real>("dFeV")),
   _dFeI(declareProperty<Real>("dFeI")),
   _dCrV_c(declareProperty<Real>("dCrV_c")),
   _dCrI_c(declareProperty<Real>("dCrI_c")),
   _dFeV_c(declareProperty<Real>("dFeV_c")),
   _dFeI_c(declareProperty<Real>("dFeI_c")),
   _kT(declareProperty<Real>("kT")),
   _Temp(declareProperty<Real>("Temp")),
   _Zg(declareProperty<Real>("Zg")),
   _dgHVaUI(declareProperty<Real>("dgHVaUI")),
   _M_Cr(declareProperty<Real>("M_Cr")),
   _M_Va(declareProperty<Real>("M_Va")),
   _M_SI(declareProperty<Real>("M_SI")),
   _M_Va_sink(declareProperty<Real>("M_Va_sink")),
   _M_SI_sink(declareProperty<Real>("M_SI_sink")),
   _M_Va_int(declareProperty<Real>("M_Va_int")),
   _M_SI_int(declareProperty<Real>("M_SI_int")),
   _DM_Cr(declareProperty<Real>("DM_Cr")),
   _DM_Va(declareProperty<Real>("DM_Va")),
   _DM_SI(declareProperty<Real>("DM_SI")),
   _Dgrad_Mnp_Cr(declareProperty<RealGradient>("Dgrad_Mnp_Cr")),
   _Dgrad_Mnp_Va(declareProperty<RealGradient>("Dgrad_Mnp_Va")),
   _Dgrad_Mnp_SI(declareProperty<RealGradient>("Dgrad_Mnp_SI")),
   _Dgrad_Mngp_Cr(declareProperty<Real>("Dgrad_Mngp_Cr")),
   _Dgrad_Mngp_Va(declareProperty<Real>("Dgrad_Mngp_Va")),
   _Dgrad_Mngp_SI(declareProperty<Real>("Dgrad_Mngp_SI")),
   _grad_M_Cr(declareProperty<RealGradient>("grad_M_Cr")),
   _grad_M_Va(declareProperty<RealGradient>("grad_M_Va")),
   _grad_M_SI(declareProperty<RealGradient>("grad_M_SI")),
   _kappa_cFe(declareProperty<Real>("kappa_cFe")),
   _kappa_cCr(declareProperty<Real>("kappa_cCr")),
   _kappa_cVa(declareProperty<Real>("kappa_cVa")),
   _kappa_cSI(declareProperty<Real>("kappa_cSI")),
   _cFe(declareProperty<Real>("cFe")),
   _grad_cFe(declareProperty<RealGradient>("grad_cFe")),
   _ICrVa(declareProperty<Real>("ICrVa")),
   _ICrSI(declareProperty<Real>("ICrSI")),
   _IVaVa(declareProperty<Real>("IVaVa")),
   _LogC(declareProperty<Real>("LogC")),
   _LogC2(declareProperty<Real>("LogC2")),
   _LogTol(declareProperty<Real>("LogTol")),
   _LNcVa(declareProperty<Real>("LNcVa")),
   _LNcSI(declareProperty<Real>("LNcSI")),
   _grad_LNcVa(declareProperty<RealGradient>("grad_LNcVa")),
   _grad_LNcSI(declareProperty<RealGradient>("grad_LNcSI"))
{
 
}

void
FeCrVaSIbulk::computeProperties()
{

  for(_qp=0; _qp < _qrule->n_points(); ++_qp)
    {
      
      _Temp[_qp] = _temp;       
      _kT[_qp] = 8.617343e-5*_temp; //Boltzmann constant in eV/K x temp
      _Zg[_qp] = _Z;
      _dgHVaUI[_qp] = _dHVaUI;
      
      _cFe[_qp]=1.0-_cCr[_qp];
      _grad_cFe[_qp]=-_grad_cCr[_qp];

      // Diffusivities assume no concentration is less than 1e-22
      Real Ca = 0.0; //
      Real Cb = 1e-6; // <- This one is in a denom...
      Real Cv = 0.0; //
      Real Ci = 0.0; //
      if(_cCr[_qp]>Ca)
      {
        Ca = _cCr[_qp];
      }
      if(_cFe[_qp]>Cb)
      {
        Cb = _cFe[_qp];
      }
      if(_cVa[_qp]>Cv)
      {
        Cv = _cVa[_qp];
      }
      if(_cSI[_qp]>Ci)
      {
        Ci = _cSI[_qp];
      }
      if(_cCr[_qp]>1.0)
      {
        Ca = 1.0;
      }
      if(_cFe[_qp]>1.0)
      {
        Cb = 1.0;
      }
      if(_cVa[_qp]>1.0)
      {
        Cv = 1.0;
      }
      if(_cSI[_qp]>1.0)
      {
        Ci = 1.0;
      }
      
      Real _LCrCr_V;
      Real _LFeFe_V;
      Real _LCrFe_V;
      Real _LCrCr_I;
      Real _LFeFe_I;
      Real _LCrFe_I;

      _LCrCr_V = _Laav1 ;
      _LFeFe_V = _Lbbv1 ;
      _LCrFe_V = _Labv1 ;
      _LCrCr_I = _Laai1 ;
      _LFeFe_I = _Lbbi1 ;
      _LCrFe_I = _Labi1 ;
      
      // ORIGINAL:
      //_dCrV[_qp] = (_LCrCr_V+_LCrFe_V) / (Ca*Cv);
      //_dCrV_c[_qp] = _LCrCr_V/(Ca*Cv) - _LCrFe_V/(Cb*Cv); // Ignore zeta term
      //_dCrI[_qp] = (_LCrCr_I+_LCrFe_I) / (Ca*Ci);
      //_dCrI_c[_qp] = _LCrCr_I/(Ca*Ci) - _LCrFe_I/(Cb*Ci); // Ignore zeta term
      //_dFeV[_qp] = (_LFeFe_V+_LCrFe_V) / (Cb*Cv);
      //_dFeV_c[_qp] = _LFeFe_V/(Cb*Cv) - _LCrFe_V/(Ca*Cv); // Ignore zeta term
      //_dFeI[_qp] = (_LFeFe_I+_LCrFe_I) / (Cb*Ci);
      //_dFeI_c[_qp] = _LFeFe_I/(Cb*Ci) - _LCrFe_I/(Ca*Ci); // Ignore zeta term
      //_M_Cr[_qp]  =  (Ca/_kT[_qp])*(_dCrV_c[_qp]*Cv+_dCrI_c[_qp]*Ci);
      //_M_Va[_qp]  =  (Ca*Cv/_kT[_qp])*(_dFeV_c[_qp]-_dCrV_c[_qp]);
      //_M_SI[_qp]  =  (Ca*Ci/_kT[_qp])*(_dCrI_c[_qp]-_dFeI_c[_qp]);
        
      // MULTIPLY L's BY DEFECT CONCENTRATION and usual d's by Ca:
      //  ..(Allows C_Cr to go to zero)
      _dCrV[_qp]   = _LCrCr_V + _LCrFe_V ;
      _dCrV_c[_qp] = _LCrCr_V - _LCrFe_V*(Ca/Cb) ; // Ignore zeta term
      _dCrI[_qp]   = _LCrCr_I + _LCrFe_I ;
      _dCrI_c[_qp] = _LCrCr_I - _LCrFe_I*(Ca/Cb) ; // Ignore zeta term
      _dFeV[_qp]   = _LFeFe_V + _LCrFe_V ;
      _dFeV_c[_qp] = _LFeFe_V*(Ca/Cb) - _LCrFe_V ; // Ignore zeta term
      _dFeI[_qp]   = _LFeFe_I + _LCrFe_I ;
      _dFeI_c[_qp] = _LFeFe_I*(Ca/Cb) - _LCrFe_I ; // Ignore zeta term
      
      // Hack TEST:
      //_dCrV_c[_qp] = _dCrV[_qp] ;
      //_dCrI_c[_qp] = _dCrI[_qp] ;
      //_dFeV_c[_qp] = _dFeV[_qp] ;
      //_dFeI_c[_qp] = _dFeI[_qp] ;
      
      _M_Cr[_qp]  =  (1.0/_kT[_qp])*(_dCrV_c[_qp]*Cv+_dCrI_c[_qp]*Ci);
      _M_Va[_qp]  =  (Cv/_kT[_qp])*(_dFeV_c[_qp]-_dCrV_c[_qp]);
      _M_SI[_qp]  =  (Ci/_kT[_qp])*(_dCrI_c[_qp]-_dFeI_c[_qp]);
      
      // These derivatives are not correct yet:
      _DM_Cr[_qp] =  0.0; //(1.0/_kT[_qp])*(_dCrV_c[_qp]*Cv+_dCrI_c[_qp]*Ci);
      _DM_Va[_qp] =  0.0; //(Ca*1.0/_kT[_qp])*(_dFeV_c[_qp]-_dCrV_c[_qp]);
      _DM_SI[_qp] =  0.0; //(Ca*1.0/_kT[_qp])*(_dCrI_c[_qp]-_dFeI_c[_qp]);
      _Dgrad_Mnp_Cr[_qp] = 0.0; //_DM_Cr[_qp] * _grad_cCr[_qp];
      _Dgrad_Mnp_Va[_qp] = 0.0; //_DM_Va[_qp] * _grad_cVa[_qp];
      _Dgrad_Mnp_SI[_qp] = 0.0; //_DM_SI[_qp] * _grad_cSI[_qp];
      _Dgrad_Mngp_Cr[_qp] = 0.0 ;
      _Dgrad_Mngp_Va[_qp] = 0.0 ;
      _Dgrad_Mngp_SI[_qp] = 0.0 ;
      _grad_M_Cr[_qp] = 0.0;
      _grad_M_Va[_qp] = 0.0;
      _grad_M_SI[_qp] = 0.0;
      
      _kappa_cFe[_qp] = _kappa_cfe;
      _kappa_cCr[_qp] = _kappa_ccr;
      _kappa_cVa[_qp] = _kappa_cva;
      _kappa_cSI[_qp] = _kappa_csI;
      
      _DrecombVa[_qp] = _LCrCr_V + _LFeFe_V + 2.0*_LCrFe_V ;
      _DrecombSI[_qp] = _LCrCr_I + _LFeFe_I + 2.0*_LCrFe_I ;
      _Drecomb[_qp] = _DrecombVa[_qp] + _DrecombSI[_qp];
      
      _M_Va_sink[_qp] = _DrecombVa[_qp]/_kT[_qp];
      _M_SI_sink[_qp] = _DrecombSI[_qp]/_kT[_qp];
      
      _M_Va_int[_qp] = (_dCrV[_qp]+_dFeV[_qp])*Cv/_kT[_qp];
      _M_SI_int[_qp] = (_dCrI[_qp]+_dFeI[_qp])*Ci/_kT[_qp];

      _ICrSI[_qp] = _ICrsI;
      _ICrVa[_qp] = _ICrva;
      _IVaVa[_qp] = _IVava;

      //Cutoffs for mathematical functions                                                                                             
      _LogC[_qp] = 1e-10;
      _LogC2[_qp] = 1e-10;
      _LogTol[_qp] = 1e-10;
      
      _LNcVa[_qp] = MathUtils::poly4Log(_cVa[_qp],_LogTol[_qp],0);
      _LNcSI[_qp] = MathUtils::poly4Log(_cSI[_qp],_LogTol[_qp],0);
 
    }
    // FOR DEBUGGING:
    //std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
    //std::cout << "_LCrCr_V = " << _LCrCr_V << std::endl;
    //std::cout << "_LCrFe_V = " << _LCrFe_V << std::endl;
    //std::cout << "_LCrCr_I = " << _LCrCr_I << std::endl;
    //std::cout << "_LCrFe_I = " << _LCrFe_I << std::endl;
    //std::cout << "_dCrV    = " << _dCrV[0] << std::endl;
    //std::cout << "_dCrI    = " << _dCrI[0] << std::endl;
    //std::cout << "_dCrV_c  = " << _dCrV_c[0] << std::endl;
    //std::cout << "_dCrI_c  = " << _dCrI_c[0] << std::endl;
    //std::cout << "_M_Cr    = " << _M_Cr[0] << std::endl;
    //std::cout << "_M_SI    = " << _M_SI[0] << std::endl;
    //std::cout << " ~~~~~~~~~~~~~~~~~~~~~~~~~~ " << std::endl;
}
