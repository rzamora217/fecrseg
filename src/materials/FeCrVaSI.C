#include "FeCrVaSI.h"
#include "MathUtils.h"
// libMesh includes
#include "libmesh/quadrature.h"
using namespace MathUtils;

template<>
InputParameters validParams<FeCrVaSI>()
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
  params.addParam<Real>("f0s", 0.125,"The GB energy constant ");
  params.addRequiredParam<Real>("M0GB","GB mobility constant in m^4/(J-s)");
  params.addRequiredParam<Real>("QGB","GB energy of migration in eV");
  params.addRequiredParam<Real>("sigma_GB","GB Surface energy in J/m^2");
  params.addRequiredCoupledVarWithAutoBuild("v", "var_name_base", "op_num", "Array of coupled variables");
  params.addRequiredParam<unsigned int>("crys_num", "number of grains");
  //params.addRequiredParam<unsigned int>("op_num", "number of grains, used to generate vector v");
  //params.addRequiredParam<std::string>("var_name_base","base for variable names");
  params.addRequiredCoupledVar("cCr","Concentration");
  params.addRequiredCoupledVar("cVa","Concentration");
  params.addRequiredCoupledVar("cSI","Concentration");
  params.addParam<Real>("tstop", 0.0005,"Stop time for grain growth ");
  params.addCoupledVar("v", "Array of coupled variables");

  return params;
}

FeCrVaSI::FeCrVaSI(const InputParameters & parameters)
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
   _f0s(getParam<Real>("f0s")),
   _M0GB(getParam<Real>("M0GB")),
   _QGB(getParam<Real>("QGB")),
   _sigma_GB(getParam<Real>("sigma_GB")),
   _tstop(getParam<Real>("tstop")),
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
   _kb(declareProperty<Real>("kb")),
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
   _kappa_GB(declareProperty<Real>("kappa_op")),
   _gamma(declareProperty<Real>("gamma_asymm")),
   _L(declareProperty<Real>("L")),
   _mu(declareProperty<Real>("mu")),
   _cFe(declareProperty<Real>("cFe")),
   _grad_cFe(declareProperty<RealGradient>("grad_cFe")),
   _ICrVa(declareProperty<Real>("ICrVa")),
   _ICrSI(declareProperty<Real>("ICrSI")),
   _IVaVa(declareProperty<Real>("IVaVa")),
   _vC_i(declareProperty<std::vector<Real> >("vC_i")),
   _vm_i(declareProperty<std::vector<Real> >("vm_i")),
   _vb_i(declareProperty<std::vector<Real> >("vb_i")),
   _vVaC_i(declareProperty<std::vector<Real> >("vVaC_i")),
   _vVam_i(declareProperty<std::vector<Real> >("vVam_i")),
   _vVab_i(declareProperty<std::vector<Real> >("vVab_i")),
   _vUIC_i(declareProperty<std::vector<Real> >("vUIC_i")),
   _vUIm_i(declareProperty<std::vector<Real> >("vUIm_i")),
   _vUIb_i(declareProperty<std::vector<Real> >("vUIb_i")),
   _LogC(declareProperty<Real>("LogC")),
   _LogC2(declareProperty<Real>("LogC2")),
   _LogTol(declareProperty<Real>("LogTol")),
   _tgrad_corr_mult(declareProperty<Real>("tgrad_corr_mult")),
   _n(coupledComponents("v"))

{

  //std::cout<<" CHECK A"<<std::endl;
  
  //_n = coupledComponents("v");
  _vals.resize(_n);
  _grad_vals.resize(_n);

  for (unsigned int i=0; i<_n; ++i)
    {
    _vals[i] = &coupledValue("v", i);
    _grad_vals[i] = &coupledGradient("v", i);
    }


  Real JtoeV = 6.24150974e18;// joule to eV conversion
  //Convert parameter units
  _sigma_GB *= JtoeV*(1e-9*1e-9);//convert of eV/nm^2
  _M0GB *= 1.0e0/(JtoeV*(1.e-9*1.e-9*1.e-9*1.e-9));//Convert to nm^4/(eV*s)
 
  //std::cout<<" CHECK B"<<std::endl;
 
}

void
FeCrVaSI::computeProperties()
{

  //std::cout<<" CHECK #1"<<std::endl;
 
  Real vC_i[] = 
    {  
       0.1584, 0.1584, 0.1584, 0.1584, 0.1584, 0.1584
    };
  Real vm_i[] =
    {
       0.04, 0.04, 0.04, 0.04, 0.04, 0.04
    };
  Real vb_i[] = 
    {
       -0.05, -0.05, -0.05, -0.05, -0.05, -0.05
    };
  Real gbw_i[] = 
    {
       4.704, 4.704, 4.704, 4.704, 4.704, 4.704
    };

//For vacancies
 Real vVaC_i[] =
    {
       1.424, 1.424, 1.424, 1.424, 1.424, 1.424
   };
  Real vVam_i[] =
    {
       1.9, 1.9, 1.9, 1.9, 1.9, 1.9
    };
  Real vVab_i[] =
    {
       -7.8, -7.8, -7.8, -7.8, -7.8, -7.8
    };

  //For inersttials                                                                                                                         
 Real vUIC_i[] =
   {
       3.012, 3.012, 3.012, 3.012, 3.012, 3.012
   };
  Real vUIm_i[] =
    {
       2.1, 2.1, 2.1, 2.1, 2.1, 2.1
    };
  Real vUIb_i[] =
    {
       -8.2, -8.2, -8.2, -8.2, -8.2, -8.2
    };

  //std::cout<<" CHECK #2"<<std::endl;

  for(_qp = 0; _qp < _qrule->n_points(); _qp++)
    {
      
      _kb[_qp] = 8.617343e-5;//Boltzmann constant in eV/K
   
      _Temp[_qp] = _temp;       
      _kT[_qp] = _kb[_qp]*_temp;
      _Zg[_qp] = _Z;
      _dgHVaUI[_qp] = _dHVaUI;
      
      _cFe[_qp]=1.0-_cCr[_qp];
      _grad_cFe[_qp]=-_grad_cCr[_qp];

      Real f_i = 0.0;
      RealGradient grad_f_i = 0.0;
      Real w_GB = 0.0;
      unsigned int ind=0;
      Real smphi = 0.0;
      
      //std::cout<<" CHECK #3"<<std::endl;
      
      for (unsigned int i=0; i<_n; ++i)
	    for (unsigned int j=i+1; j<_n; ++j)
	    {
	      f_i += -8.0*(*_vals[i])[_qp]*(*_vals[i])[_qp]*(*_vals[j])[_qp]*(*_vals[j])[_qp];
	      grad_f_i += -8.0*2.0*(*_vals[i])[_qp]*(*_grad_vals[i])[_qp]*(*_vals[j])[_qp]*(*_vals[j])[_qp];
	      grad_f_i += -8.0*(*_vals[i])[_qp]*(*_vals[i])[_qp]*2.0*(*_vals[j])[_qp]*(*_grad_vals[j])[_qp];
	      
	      w_GB += gbw_i[ind]*(*_vals[i])[_qp]*(*_vals[i])[_qp]*(*_vals[j])[_qp]*(*_vals[j])[_qp];
	      smphi += (*_vals[i])[_qp]*(*_vals[i])[_qp]*(*_vals[j])[_qp]*(*_vals[j])[_qp];
	  
	      ind++;
	    }
	  if (smphi > 0)
	    w_GB/=smphi;
	  else
	    w_GB = gbw_i[0];

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
      _M_Cr[_qp]  =  (1.0/_kT[_qp])*(_dCrV_c[_qp]*Cv+_dCrI_c[_qp]*Ci);
      _M_Va[_qp]  =  (Cv/_kT[_qp])*(_dFeV_c[_qp]-_dCrV_c[_qp]);
      _M_SI[_qp]  =  (Ci/_kT[_qp])*(_dCrI_c[_qp]-_dFeI_c[_qp]);
      
      // These derivatives are not correct yet:
      _DM_Cr[_qp] =  0.0 ;
      _DM_Va[_qp] =  0.0 ;
      _DM_SI[_qp] =  0.0 ;
      _Dgrad_Mnp_Cr[_qp] = 0.0 ;
      _Dgrad_Mnp_Va[_qp] = 0.0 ;
      _Dgrad_Mnp_SI[_qp] = 0.0 ;
      _Dgrad_Mngp_Cr[_qp] = 0.0 ;
      _Dgrad_Mngp_Va[_qp] = 0.0 ;
      _Dgrad_Mngp_SI[_qp] = 0.0 ;
      _grad_M_Cr[_qp] = 0.0 ;
      _grad_M_Va[_qp] = 0.0 ;
      _grad_M_SI[_qp] = 0.0 ;

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

      _vC_i[_qp].resize(6);
      _vm_i[_qp].resize(6);
      _vb_i[_qp].resize(6); 
      
      _vVaC_i[_qp].resize(6);
      _vVam_i[_qp].resize(6);
      _vVab_i[_qp].resize(6);
      
      _vUIC_i[_qp].resize(6);
      _vUIm_i[_qp].resize(6);
      _vUIb_i[_qp].resize(6);

      for (unsigned int i=0; i<6; ++i)
	{
	  _vC_i[_qp][i] = vC_i[i];
	  _vm_i[_qp][i] = vm_i[i]; 
	  _vb_i[_qp][i] = vb_i[i]; 
	}
      
      for (unsigned int i=0; i<6; ++i)
     	{
	  _vVaC_i[_qp][i] = vVaC_i[i];
      	  _vVam_i[_qp][i] = vVam_i[i];
      	  _vVab_i[_qp][i] = vVab_i[i];
	}

      for (unsigned int i=0; i<6; ++i)
        {
          _vUIC_i[_qp][i] = vUIC_i[i];
          _vUIm_i[_qp][i] = vUIm_i[i];
          _vUIb_i[_qp][i] = vUIb_i[i];
        }

      //For GB
       Real M_GB=0.0;
       if (_t>_tstop) {
         M_GB=0.0;
        }
       else {
        M_GB = _M0GB*std::exp(-_QGB/(_kb[_qp]*_temp));
        }
      _L[_qp] = 4.0/3.0*M_GB/w_GB;
      _kappa_GB[_qp] = 3.0/4.0*_sigma_GB*w_GB;
      _mu[_qp] = 3.0/4.0*1/_f0s*_sigma_GB/w_GB;
      _tgrad_corr_mult[_qp] = _mu[_qp]*9.0/8.0;
      _gamma[_qp] = 1.5;

      //Cutoffs for mathematical functions                                                                                             
      _LogC[_qp] = 1e-10;
      _LogC2[_qp] = 1e-10;
      _LogTol[_qp] = 1e-10;
 
    }
}
