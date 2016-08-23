#!/usr/bin/env python

import os
import string
import time
import random
import numpy

def writeinput(casename,temp,conc,xdim,nx,LCrCr_V,LFeFe_V,LCrFe_V,LCrCr_I,LFeFe_I,LCrFe_I,nsteps,dt,GBbc,kappa_ccr,kappa_cva,kappa_csI,Xv,Xi,use_split,dose_rate,end_time,use_dirichlet):

  generatedmesh = True

  # Write New Input File:
  filename = casename+'_'+GBbc+'.i'
  f_out = open(filename, 'w')
  f_out.write('# MARMOT INPUT FILE: %s \n' %(filename))
  f_out.write('\n')
  f_out.write('[Mesh]\n')
  
  if generatedmesh:
    f_out.write('  type = GeneratedMesh\n')
    f_out.write('  dim = 2\n')
    if GBbc == 'ideal':
      f_out.write('  nx = %d\n' %(nx))
    else:
      f_out.write('  nx = %d\n' %(2*nx))
    f_out.write('  ny = 1\n')
    f_out.write('  nz = 0\n')
    f_out.write('  xmin = 0\n')
    if GBbc == 'ideal':
      f_out.write('  xmax = %d\n' %(xdim))
    else:
      f_out.write('  xmax = %d\n' %(2*xdim))
    f_out.write('  ymin = 0\n')
    f_out.write('  ymax = 1\n')
    f_out.write('  zmin = 0\n')
    f_out.write('  zmax = 0\n')
    f_out.write('  elem_type = QUAD4\n')
  else:
    f_out.write('  type = FileMesh\n')
    f_out.write('  file = ricksmesh.e\n')
    f_out.write('  block_id = \'9\'\n')
    f_out.write('  boundary_id = \'1 2\'\n')
    f_out.write('  boundary_name = \'left right\'\n')
  f_out.write('[]\n')
  f_out.write('\n')
  if GBbc != 'ideal':
    f_out.write('[GlobalParams]\n')
    f_out.write('  crys_num = 2\n')
    f_out.write('  op_num = 2\n')
    f_out.write('  var_name_base = gr\n')
    f_out.write('[]\n')
    f_out.write('\n')
  f_out.write('[Variables]\n')
  f_out.write('  [./cCr]\n')
  f_out.write('    order = THIRD\n')
  f_out.write('    family = HERMITE\n')
  f_out.write('    scaling=1.00E-00\n')
  f_out.write('    block = 0\n')
  f_out.write('  [../]\n')
  if use_split:
    f_out.write('  [./wCr]\n')
    f_out.write('    order = THIRD #FIRST\n')
    f_out.write('    family = HERMITE #LAGRANGE\n')
    f_out.write('    scaling=1.00E-00\n')
    f_out.write('    block = 0\n')
    f_out.write('  [../]\n')
  f_out.write('  [./cVa]\n')
  f_out.write('    order = THIRD #FIRST\n')
  f_out.write('    family = HERMITE #LAGRANGE\n')
  f_out.write('    scaling=1.00E-00\n')
  f_out.write('    block = 0\n')
  f_out.write('  [../]\n')
  f_out.write('  [./cSI]\n')
  f_out.write('    order = THIRD #FIRST\n')
  f_out.write('    family = HERMITE #LAGRANGE\n')
  f_out.write('    scaling=1.00E-00\n')
  f_out.write('    block = 0\n')
  f_out.write('  [../]\n')
  if GBbc != 'ideal':
    f_out.write('  [./PolycrystalVariables]\n')
    f_out.write('  [../]\n')
    f_out.write('\n')
  
  f_out.write('[]\n')
  f_out.write('\n')
  f_out.write('[Materials]\n')
  f_out.write('  active = \'LANL free_energy\'\n')
  f_out.write('  [./LANL]\n')
  f_out.write('\n')
  if GBbc == 'ideal':
    f_out.write('    type = FeCrVaSIbulk\n')
  else:
    f_out.write('    type = FeCrVaSI\n')
  f_out.write('\n')
  f_out.write('    cCr = cCr\n')
  f_out.write('    cVa = cVa\n')
  f_out.write('    cSI = cSI\n')
  f_out.write('\n')
  f_out.write('    kappa_cfe = 1.0\n')
  f_out.write('    kappa_cva = %s\n' %(kappa_cva))
  f_out.write('    kappa_ccr = %s\n' %(kappa_ccr))
  f_out.write('    kappa_csI = %s\n' %(kappa_csI))
  f_out.write('\n')
  if GBbc != 'ideal':
    f_out.write('    M0GB = 1.0e-4 # m^4/(J*s); Increase for smalle time steps 2.5e-6\n')
    f_out.write('    QGB = 1.3 #65 #2.0 #1.6 # eV\n')
    f_out.write('    sigma_GB = 0.708 # J/m^2\n')
    # UO2 (2000K) -> M0GB = 2.5e-6 # m^4/(J*s); Increase for smalle time steps 2.5e-6
    # UO2 (2000K) -> QGB = 2.3 # eV
    f_out.write('    crys_num = 2\n')
    f_out.write('    op_num = 2\n')
    f_out.write('    var_name_base = gr\n')
    f_out.write('    tstop = 1.0e-0\n')
    f_out.write('\n')
  f_out.write('    ICrva = 0.0 #-0.048 # eV\n')
  f_out.write('    ICrsI = 0.0 #-0.065 # eV\n')
  f_out.write('    IVava = 0.0 #-0.049 # eV\n')
  f_out.write('\n')
  f_out.write('    temp = %s # K\n' %(temp))
  f_out.write('    Laav1 = %s\n' %(LCrCr_V))
  f_out.write('    Lbbv1 = %s\n' %(LFeFe_V))
  f_out.write('    Labv1 = %s\n' %(LCrFe_V))
  f_out.write('    Laai1 = %s\n' %(LCrCr_I))
  f_out.write('    Lbbi1 = %s\n' %(LFeFe_I))
  f_out.write('    Labi1 = %s\n' %(LCrFe_I))
  #f_out.write('    drecombVa = %s\n' %(drecombVa))
  #f_out.write('    drecombSI = %s\n' %(drecombSI))
  f_out.write('\n')
  f_out.write('    Z = 1.0\n')
  f_out.write('    dHVaUI = 1.2 #eV\n')
  f_out.write('    block = 0\n')
  f_out.write('  [../]\n')
  f_out.write('  [./free_energy]\n')
  f_out.write('    type = DerivativeParsedMaterial\n')
  f_out.write('    block = 0\n')
  f_out.write('    f_name = F\n')
  f_out.write('    args = \'cCr cVa cSI\'\n')
  f_out.write('    constant_names       = \'T0  y   kB\'\n')
  f_out.write('    constant_expressions = \'410 %s 8.6173324e-5\'\n' %(temp))
  #f_out.write('    third_derivatives = false\n')
  f_out.write('    derivative_order = 2\n')
  f_out.write('    enable_jit = true\n')
  f_out.write('\n')
  f_out.write('    # G(CFe) from Mathematica\n')
  if False:
    f_out.write('    function = \'y*kB*cCr*plog(cCr,0.0001)\'\n')
  else:
    f_out.write('    function = \'(cCr*-4.16723234898182 + (1 - cCr)*-3.85953051397515 + cCr*(1 - cCr)*(0.380577499303532 + 0.0885190880468175*(1 - 2*cCr) + -0.0549325534161349*(1 - 2*cCr)^2 + 0.203752174307515*(1 - 2*cCr)^3 + -0.164028794268281*(1 - 2*cCr)^4 + -0.0172550241855144*(1 - 2*cCr)^5)) * y/T0\n')
    f_out.write('    + (cCr*-4.12297413272388 + (1 - cCr)*-3.83660056975999 + cCr*(1 - cCr)*(0.385512095260579 + 0.0962430189496054*(1 - 2*cCr) + -0.041249704177587*(1 - 2*cCr)^2 + 0.194439246552959*(1 - 2*cCr)^3 + -0.195412847295217*(1 - 2*cCr)^4 + 0.00967038578662529*(1 - 2*cCr)^5))\n')
    f_out.write('    -(0.000263056717498972 + 4.614980081531*10^-5*cCr + -4.75235048526914*10^-5*cCr^2 + 9.35759929354588*10^-6*cCr^3)*y*log(y)\n')
    f_out.write('    - (3.04676203180853*10^-9 + -2.07225774483557*10^-8*cCr + 3.55582178830517*10^-8*cCr^2 + -2.70425743485173*10^-8*cCr^3)*y^2\n')
    f_out.write('    + (-1.51823088659839*10^-13 + 5.18553402098699*10^-12*cCr + -4.56309143596694*10^-12*cCr^2 + 1.08597105154957*10^-11*cCr^3)/2*y^3\n')
    f_out.write('    + (-(cCr*-4.12297413272388 + (1 - cCr)*-3.83660056975999 +\n')
    f_out.write('     cCr*(1 - cCr)*(0.385512095260579 +\n')
    f_out.write('         0.0962430189496054*(1 - 2*cCr) + -0.041249704177587*(1 - 2*cCr)^2 +\n')
    f_out.write('             0.194439246552959*(1 - 2*cCr)^3 + -0.195412847295217*(1 - 2*cCr)^4 +\n')
    f_out.write('                 0.00967038578662529*(1 - 2*cCr)^5))/T0 + (3.04676203180853*10^-9 + -2.07225774483557*10^-8*cCr +\n')
    f_out.write('                  3.55582178830517*10^-8*cCr^2 + -2.70425743485173*10^-8*cCr^3)*T0 + (0.000263056717498972 +\n')
    f_out.write('                   4.614980081531*10^-5*cCr + -4.75235048526914*10^-5*cCr^2 +\n')
    f_out.write('                    9.35759929354588*10^-6*cCr^3)*log(T0) + (-1.51823088659839*10^-13 +\n')
    f_out.write('                     5.18553402098699*10^-12*cCr + -4.56309143596694*10^-12*cCr^2 +\n')
    f_out.write('                      1.08597105154957*10^-11*cCr^3)/2*T0^2)*y + y*kB*(cCr*plog(cCr,0.0001) + (1 - cCr)*plog(1 - cCr,0.0001))\'\n')
  f_out.write('  [../]\n')
  f_out.write('[]\n')
  f_out.write('\n')
  f_out.write('[Kernels]\n')
  
  if GBbc == 'ideal':
    fcut = 0.1
    if use_split:
      f_out.write('  active = \'cCr_CHRes wCr_CHRes time_cCr cVa_Diffusion ie_cVa cSI_Diffusion UISource VaSource Rate_FeCrSI Rate_FeCrVa\'\n')
    else:
      #f_out.write('  active = \'cCr_CHParsed cCr_CHint cCr_DiffusionVa cCr_DiffusionSI ie_cCr cVa_CHParsed cVa_CHint cVa_Diffusion ie_cVa cSI_CHParsed cSI_CHint cSI_Diffusion ie_cSI UISource VaSource Rate_FeCrSI Rate_FeCrVa\'\n')
      #f_out.write('  active = \'cCr_CHParsed cCr_CHint cCr_DiffusionVa cCr_DiffusionSI ie_cCr cVa_CHParsed cVa_Diffusion ie_cVa cSI_CHParsed cSI_Diffusion UISource VaSource Rate_FeCrSI Rate_FeCrVa\'\n')
      #f_out.write('  active = \'cCr_CHParsed cCr_CHint cCr_DiffusionVa cCr_DiffusionSI ie_cCr cVa_Diffusion ie_cVa cSI_Diffusion UISource VaSource Rate_FeCrSI Rate_FeCrVa\'\n')
      #f_out.write('  active = \'ie_cCr ie_cVa ie_cSI UISource VaSource Rate_FeCrSI Rate_FeCrVa\'\n')
    
      #f_out.write('  active = \'ie_cCr cVa_Diffusion ie_cVa ie_cSI VaSource\'\n')
    
      #f_out.write('  active = \'ie_cCr cVa_Diffusion ie_cVa cSI_Diffusion ie_cSI UISource VaSource Rate_FeCrSI Rate_FeCrVa\'\n')

      #f_out.write('  active = \'cCr_CHParsed cCr_CHint cCr_DiffusionVa cCr_DiffusionSI ie_cCr cVa_Diffusion ie_cVa cSI_Diffusion ie_cSI UISource VaSource Rate_FeCrSI Rate_FeCrVa\'\n')
      
      #f_out.write('  active = \'cCr_CHParsed cCr_CHint cCr_DiffusionVa cCr_DiffusionSI ie_cCr cVa_CHParsed cVa_Diffusion ie_cVa cSI_CHParsed cSI_Diffusion ie_cSI UISource VaSource Rate_FeCrSI Rate_FeCrVa\'\n')
      if use_dirichlet:
        f_out.write('  active = \'cCr_CHParsed cCr_CHint cCr_DiffusionVa cCr_DiffusionSI ie_cCr cVa_CHParsed cVa_Diffusion cVa_CHint ie_cVa cSI_CHParsed cSI_Diffusion cSI_CHint ie_cSI UISource VaSource Rate_FeCrSI Rate_FeCrVa\'\n')
      else:
        #f_out.write('  active = \'cCr_CHParsed cCr_CHint cCr_DiffusionVa cCr_DiffusionSI ie_cCr cVa_CHParsed cVa_Diffusion cVa_CHint ie_cVa cSI_CHParsed cSI_Diffusion cSI_CHint ie_cSI UISource VaSource Rate_FeCrSI Rate_FeCrVa Va_idealsink SI_idealsink\'\n')
    
        f_out.write('  active = \'cCr_CHParsed cCr_CHint cCr_DiffusionVa cCr_DiffusionSI ie_cCr cVa_CHParsed cVa_Diffusion cVa_CHint ie_cVa cSI_CHParsed cSI_Diffusion cSI_CHint ie_cSI UISource VaSource Rate_FeCrSI Rate_FeCrVa Va_GBsink SI_GBsink\'\n')

    f_out.write('  [./Va_idealsink]\n')
    f_out.write('    type = FeCrIdealSink\n')
    f_out.write('    variable = cVa\n')
    f_out.write('    c = cSI\n')
    f_out.write('    xdim = %f\n' %(xdim))
    f_out.write('    dose_rate = %f\n' %(dose_rate))
    f_out.write('    fcut = %f\n' %(fcut))
    f_out.write('    diffusivitySink = DrecombVa\n')
    f_out.write('  [../]\n')
    f_out.write('  [./SI_idealsink]\n')
    f_out.write('    type = FeCrIdealSink\n')
    f_out.write('    variable = cSI\n')
    f_out.write('    c = cVa\n')
    f_out.write('    xdim = %f\n' %(xdim))
    f_out.write('    dose_rate = %f\n' %(dose_rate))
    f_out.write('    fcut = %f\n' %(fcut))
    f_out.write('    diffusivitySink = DrecombSI\n')
    f_out.write('  [../]\n')
    f_out.write('  [./Cr_GBsink]\n')
    f_out.write('    type = CHRadSinkFeCrNoPF\n')
    f_out.write('    xdim = %f\n' %(xdim))
    f_out.write('    variable = cCr\n')
    f_out.write('    mob_name = M_Cr\n')
    f_out.write('    LogC_name = LogC\n')
    f_out.write('    LogTol_name = LogTol\n')
    f_out.write('  [../]\n')
    f_out.write('  [./Va_GBsink]\n')
    f_out.write('    type = CHRadSinkFeVaNoPF\n')
    f_out.write('    xdim = %f\n' %(xdim))
    f_out.write('    variable = cVa\n')
    f_out.write('    mob_name = M_Va_sink\n')
    f_out.write('    LogC_name = LogC\n')
    f_out.write('    LogTol_name = LogTol\n')
    f_out.write('  [../]\n')
    f_out.write('  [./SI_GBsink]\n')
    f_out.write('    type = CHRadSinkFeSINoPF\n')
    f_out.write('    xdim = %f\n' %(xdim))
    f_out.write('    variable = cSI\n')
    f_out.write('    mob_name = M_SI_sink\n')
    f_out.write('    LogC_name = LogC\n')
    f_out.write('    LogTol_name = LogTol\n')
    f_out.write('  [../]\n')
  else:
    fcut = 0.0
    #f_out.write('  active = \'PolycrystalKernel cCr_CHParsed cCr_CHint cCr_DiffusionVa cCr_DiffusionSI ie_cCr cVa_CHParsed cVa_Diffusion ie_cVa cSI_CHParsed cSI_Diffusion ie_cSI UISource VaSource Rate_FeCrSI Rate_FeCrVa cui_void2 cui_void3\'\n')
    
    #f_out.write('  active = \'PolycrystalKernel cCr_CHParsed cCr_CHint cCr_DiffusionVa cCr_DiffusionSI ie_cCr cVa_CHParsed cVa_Diffusion cVa_CHint ie_cVa cSI_CHParsed cSI_Diffusion cSI_CHint ie_cSI UISource VaSource Rate_FeCrSI Rate_FeCrVa cui_void1 cui_void2 cui_void3\'\n')
    
    
    f_out.write('  active = \'PolycrystalKernel ie_cCr cVa_Diffusion ie_cVa cSI_Diffusion ie_cSI UISource VaSource Rate_FeCrSI Rate_FeCrVa cui_void2 cui_void3\'\n')
    
    f_out.write('\n')
    f_out.write('  [./PolycrystalKernel]\n')
    f_out.write('  [../]\n')
    f_out.write('\n')
    f_out.write('  [./cui_void1]\n')
    f_out.write('    type = CHRadSinkFeCr\n')
    f_out.write('    variable = cCr\n')
    f_out.write('    mob_name = M_Cr\n')
    f_out.write('    LogC_name = LogC\n')
    f_out.write('    LogTol_name = LogTol\n')
    f_out.write('  [../]\n')
    f_out.write('  [./cui_void2]\n')
    f_out.write('    type = CHRadSinkFeVa\n')
    f_out.write('    variable = cVa\n')
    f_out.write('    mob_name = M_Va_sink\n')
    f_out.write('    LogC_name = LogC\n')
    f_out.write('    LogTol_name = LogTol\n')
    f_out.write('  [../]\n')
    f_out.write('  [./cui_void3]\n')
    f_out.write('    type = CHRadSinkFeSI\n')
    f_out.write('    variable = cSI\n')
    f_out.write('    mob_name = M_SI_sink\n')
    f_out.write('    LogC_name = LogC\n')
    f_out.write('    LogTol_name = LogTol\n')
    f_out.write('  [../]\n')
    f_out.write('  [./Va_idealsink]\n')
    f_out.write('    type = FeCrIdealSink\n')
    f_out.write('    variable = cVa\n')
    f_out.write('    c = cSI\n')
    f_out.write('    xdim = %f\n' %(xdim))
    f_out.write('    dose_rate = %f\n' %(dose_rate))
    f_out.write('    fcut = %f\n' %(fcut))
    f_out.write('    diffusivitySink = DrecombVa\n')
    f_out.write('  [../]\n')
    f_out.write('  [./SI_idealsink]\n')
    f_out.write('    type = FeCrIdealSink\n')
    f_out.write('    variable = cSI\n')
    f_out.write('    c = cVa\n')
    f_out.write('    xdim = %f\n' %(xdim))
    f_out.write('    dose_rate = %f\n' %(dose_rate))
    f_out.write('    fcut = %f\n' %(fcut))
    f_out.write('    diffusivitySink = DrecombSI\n')
    f_out.write('  [../]\n')
    f_out.write('\n')
    
  f_out.write('  [./cCr_CHRes]\n')
  f_out.write('    type = SplitCHParsed\n')
  f_out.write('    variable = cCr\n')
  f_out.write('    f_name = F\n')
  f_out.write('    kappa_name = kappa_cCr\n')
  f_out.write('    w = wCr\n')
  f_out.write('  [../]\n')    
  f_out.write('  [./wCr_CHRes]\n')
  f_out.write('    type = SplitCHWRes\n')
  f_out.write('    variable = wCr\n')
  f_out.write('    mob_name = M_Cr\n')
  f_out.write('  [../]\n')
  f_out.write('  [./time_cCr]\n')
  f_out.write('    type = CoupledTimeDerivative\n')
  f_out.write('    variable = wCr\n')
  f_out.write('    v = cCr\n')
  f_out.write('  [../]\n')

  f_out.write('\n')
  f_out.write('  [./cCr_CHParsed]\n')
  #f_out.write('    type = DiffusionCHParsedCr\n')
  f_out.write('    type = CahnHilliard\n')
  f_out.write('    variable = cCr\n')
  f_out.write('    f_name = F\n')
  f_out.write('    #args = cCr\n')
  f_out.write('    mob_name = M_Cr\n')
  f_out.write('  [../]\n')
  f_out.write('  [./cCr_CHint]\n')
  #f_out.write('    type = FeCrVaSICHInt\n')
  f_out.write('    type = CHInterface\n')
  f_out.write('    variable = cCr\n')
  f_out.write('    mob_name = M_Cr\n')
  f_out.write('    kappa_name = kappa_cCr\n')
  f_out.write('    #grad_mob_name = grad_M_Cr\n')
  f_out.write('  [../]\n')
  f_out.write('  [./cCr_DiffusionVa]\n')
  f_out.write('    type = DiffusionAdd\n')
  f_out.write('    variable = cCr\n')
  f_out.write('    diffusivity = dCrV\n')
  f_out.write('    c = cVa\n')
  f_out.write('    dirscale = -1.0\n')
  f_out.write('  [../]\n')
  f_out.write('  [./cCr_DiffusionSI]\n')
  f_out.write('    type = DiffusionAdd\n')
  f_out.write('    variable = cCr\n')
  f_out.write('    diffusivity = dCrI\n')
  f_out.write('    c = cSI\n')
  f_out.write('    dirscale = 1.0\n')
  f_out.write('  [../]\n')
  f_out.write('  [./ie_cCr]\n')
  f_out.write('    type = TimeDerivative\n')
  f_out.write('    variable = cCr\n')
  f_out.write('  [../]\n')
  f_out.write('\n')
  f_out.write('  [./cVa_CHParsed]\n')
  #f_out.write('    type = DiffusionCHParsedD\n')
  f_out.write('    type = CahnHilliardD\n')
  f_out.write('    variable = cVa\n')
  f_out.write('    f_name = F\n')
  f_out.write('    args = cCr\n')
  f_out.write('    mob_name = M_Va\n')
  f_out.write('  [../]\n')
  f_out.write('  [./cVa_CHint]\n')
  f_out.write('    type = CHInterface\n')
  f_out.write('    variable = cVa\n')
  f_out.write('    mob_name = M_Va_int\n')
  f_out.write('    kappa_name = kappa_cVa\n')
  f_out.write('  [../]\n')
  f_out.write('  [./cVa_Diffusion]\n')
  f_out.write('    type = DiffusionAdd2\n')
  f_out.write('    variable = cVa\n')
  f_out.write('    diffusivityCr = dCrV\n')
  f_out.write('    diffusivityFe = dFeV\n')
  f_out.write('    cCr = cCr\n')
  f_out.write('    cVa = cVa\n')
  f_out.write('    xdim = %f\n' %(xdim))
  f_out.write('    fcut = %f\n' %(fcut))
  f_out.write('  [../]\n')
  f_out.write('  [./ie_cVa]\n')
  f_out.write('    type = TimeDerivative\n')
  f_out.write('    variable = cVa\n')
  f_out.write('  [../]\n')
  f_out.write('\n')
  f_out.write('  [./cSI_CHParsed]\n')
  #f_out.write('    type = DiffusionCHParsedD\n')
  f_out.write('    type = CahnHilliardD\n')
  f_out.write('    variable = cSI\n')
  f_out.write('    f_name = F\n')
  f_out.write('    args = cCr\n')
  f_out.write('    mob_name = M_SI\n')
  f_out.write('  [../]\n')
  f_out.write('  [./cSI_CHint]\n')
  f_out.write('    type = CHInterface\n')
  f_out.write('    variable = cSI\n')
  f_out.write('    mob_name = M_SI_int\n')
  f_out.write('    kappa_name = kappa_cCr\n')
  f_out.write('    #grad_mob_name = grad_M_SI\n')
  f_out.write('  [../]\n')
  f_out.write('  [./cSI_Diffusion]\n')
  f_out.write('    type = DiffusionAdd2\n')
  f_out.write('    variable = cSI\n')
  f_out.write('    diffusivityCr = dCrI\n')
  f_out.write('    diffusivityFe = dFeI\n')
  f_out.write('    cCr = cCr\n')
  f_out.write('    cVa = cVa\n')
  f_out.write('    xdim = %f\n' %(xdim))
  f_out.write('    fcut = %f\n' %(fcut))
  f_out.write('  [../]\n')
  f_out.write('  [./ie_cSI]\n')
  f_out.write('    type = TimeDerivative\n')
  f_out.write('    variable = cSI\n')
  f_out.write('  [../]\n')
  f_out.write('\n')
  f_out.write('  [./Rate_FeCrSI]\n')
  f_out.write('    type = FeCrDRecombRate\n')
  f_out.write('    variable = cSI\n')
  f_out.write('    diffusivityCrI = dCrI\n')
  f_out.write('    diffusivityFeI = dFeI\n')
  f_out.write('    diffusivityCrV = dCrV\n')
  f_out.write('    diffusivityFeV = dFeV\n')
  f_out.write('    cCr = cCr\n')
  f_out.write('    c = cVa\n')
  f_out.write('    cVa = cVa\n')
  f_out.write('    cSI = cSI\n')
  f_out.write('    LogC_name = LogC\n')
  f_out.write('    LogTol_name = LogTol\n')
  f_out.write('  [../]\n')
  f_out.write('  [./Rate_FeCrVa]\n')
  f_out.write('   type = FeCrDRecombRate\n')
  f_out.write('    variable = cVa\n')
  f_out.write('    diffusivityCrI = dCrI\n')
  f_out.write('    diffusivityFeI = dFeI\n')
  f_out.write('    diffusivityCrV = dCrV\n')
  f_out.write('    diffusivityFeV = dFeV\n')
  f_out.write('    cCr = cCr\n')
  f_out.write('    c = cSI\n')
  f_out.write('    cVa = cVa\n')
  f_out.write('    cSI = cSI\n')
  f_out.write('    LogC_name = LogC\n')
  f_out.write('    LogTol_name = LogTol\n')
  f_out.write('  [../]\n')
  f_out.write('\n')
  f_out.write('  [./UISource]\n')
  f_out.write('    type = UISource\n')
  f_out.write('    variable = cSI\n')
  f_out.write('    R = %f\n' %(dose_rate))
  f_out.write('    Factor = 1.0\n')
  f_out.write('  [../]\n')
  f_out.write('  [./VaSource]\n')
  f_out.write('    type = VaSource\n')
  f_out.write('    variable = cVa\n')
  f_out.write('    R = %f\n' %(dose_rate))
  f_out.write('    Factor = 1.0\n')
  f_out.write('  [../]\n')
  f_out.write('[]\n')
  f_out.write('\n')
  f_out.write('[BCs]\n')
  if GBbc == 'ideal':
    if use_dirichlet:
      f_out.write('  active = \'cVa_gb cSI_gb Periodic\'\n')
    else:
      f_out.write('  active = \'Periodic\'\n')
  else:
    f_out.write('  active = \'Periodic\'\n')
  f_out.write('  [./cVa_gb]\n')
  f_out.write('    type = DirichletBC\n')
  f_out.write('    variable = cVa\n')
  f_out.write('    boundary = \'left right\'\n')
  f_out.write('    value = %.6e\n' %(Xv))
  f_out.write('  [../]\n')
  f_out.write('  [./cVa_gbn]\n')
  f_out.write('    type = NeumannBC\n')
  f_out.write('    variable = cVa\n')
  f_out.write('    boundary = \'left right\'\n')
  f_out.write('    value = 0.0\n')
  f_out.write('  [../]\n')
  f_out.write('  [./cSI_gb]\n')
  f_out.write('    type = DirichletBC\n')
  f_out.write('    variable = cSI\n')
  f_out.write('    boundary = \'left right\'\n')
  f_out.write('    value = %.6e\n' %(Xi))
  f_out.write('  [../]\n')
  f_out.write('  [./cSI_gbn]\n')
  f_out.write('    type = NeumannBC\n')
  f_out.write('    variable = cSI\n')
  f_out.write('    boundary = \'left right\'\n')
  f_out.write('    value = 0.0\n')
  f_out.write('  [../]\n')
  f_out.write('  [./cCr_gbl]\n')
  f_out.write('    type = NeumannBC\n')
  #f_out.write('    type = DirichletBC\n')
  f_out.write('    variable = cCr\n')
  f_out.write('    boundary = \'left\'\n')
  f_out.write('    value = 0.0\n')
  f_out.write('  [../]\n')
  f_out.write('  [./cCr_gbr]\n')
  f_out.write('    type = NeumannBC\n')
  #f_out.write('    type = DirichletBC\n')
  f_out.write('    variable = cCr\n')
  f_out.write('    boundary = \'right\'\n')
  f_out.write('    value = 0.0\n')
  f_out.write('  [../]\n')
  f_out.write('  [./Periodic]\n')
  
  if GBbc == 'ideal':
    f_out.write('    [./cCr]\n')
    f_out.write('      variable = cCr\n')
    f_out.write('      auto_direction = \'x y\'\n')
    f_out.write('    [../]\n')
    f_out.write('    [./cVa]\n')
    f_out.write('      variable = cVa\n')
    f_out.write('      auto_direction = \'x y\'\n')
    f_out.write('    [../]\n')
    f_out.write('    [./cSI]\n')
    f_out.write('      variable = cSI\n')
    f_out.write('      auto_direction = \'x y\'\n')
    f_out.write('    [../]\n')
    f_out.write('  [../]\n')
    f_out.write('[]\n')
    f_out.write('\n')
  else:
    f_out.write('    [./top_bottom]\n')
    f_out.write('      primary = 0\n') #0
    f_out.write('      secondary = 2\n') #2
    f_out.write('      translation = \'0 1.0 0\'\n')
    f_out.write('    [../]\n')
    f_out.write('    [./left_right]\n')
    f_out.write('      primary = 1\n') #1
    f_out.write('      secondary = 3\n') #3
    f_out.write('      translation = \'-%d.0 0 0\'\n' %(2*xdim))
    f_out.write('    [../]\n')
    f_out.write('  [../]\n')
    f_out.write('[]\n')
    f_out.write('\n')

  f_out.write('\n')
  f_out.write('[ICs]\n')
  if GBbc == 'ideal':
    f_out.write('  active = \'cCr cVa cSI\'\n')
  else:
    f_out.write('  active = \'cCr cVa cSI PolycrystalICs\'\n')
  f_out.write('  [./cCr]\n')
  f_out.write('    type = ConstantIC\n')
  f_out.write('    variable = cCr\n')
  f_out.write('    value = %s\n' %(conc))
  f_out.write('  [../]\n')
  f_out.write('  [./cVa]\n')
  f_out.write('    type = ConstantIC\n')
  f_out.write('    variable = cVa\n')
  f_out.write('    value = %.6e\n' %(Xv))
  f_out.write('  [../]\n')
  f_out.write('  [./cSI]\n')
  f_out.write('    type = ConstantIC\n')
  f_out.write('    variable = cSI\n')
  f_out.write('    value = %.6e\n' %(Xi))
  f_out.write('  [../]\n')
  
  f_out.write('  [./cCr_Random]\n')
  f_out.write('    type = RandomIC\n')
  f_out.write('    variable = cCr\n')
  f_out.write('    min = 0.05\n')
  f_out.write('    max = 0.15\n')
  f_out.write('    block = 0\n')
  f_out.write('  [../]\n')
  f_out.write('  [./cCr_Point]\n')
  f_out.write('    type = BoundingBoxIC\n')
  f_out.write('    variable = cCr\n')
  f_out.write('    x1 = -0.001 # 4.999 #-0.001 # x coordinate of the lower left-hand corner of the box\n')
  f_out.write('    y1 = -0.001 #-0.001 #-0.001 # y coordinate of the lower left-hand corner of the box\n')
  f_out.write('    x2 =  0.001 # 5.001 # 0.001 # x coordinate of the upper right-hand corner of the box\n')
  f_out.write('    y2 =  1.001 # 1.001 # 1.001 # y coordinate of the upper right-hand corner of the box\n')
  f_out.write('    inside = 0.1\n' )
  f_out.write('    outside = %s\n' %(conc))
  f_out.write('    block = 0\n')
  f_out.write('  [../]\n')

  if GBbc != 'ideal':
    f_out.write('  [./PolycrystalICs]\n')
    f_out.write('    [./BicrystalBoundingBoxIC]\n')
    f_out.write('      x1 = 0.0\n')
    f_out.write('      y1 = 0.0\n')
    f_out.write('      x2 = %d.0\n' %(xdim))
    f_out.write('      y2 = 1.0\n')
    f_out.write('    [../]\n')
    f_out.write('  [../]\n')

  f_out.write('[]\n')
  f_out.write('\n')
  f_out.write('[Postprocessors]\n')
  f_out.write('  [./dt]\n')
  f_out.write('    type = TimestepSize\n')
  f_out.write('  [../]\n')
  f_out.write('  [./max_cCr]\n')
  f_out.write('    type = NodalMaxValue\n')
  f_out.write('    variable = cCr\n')
  f_out.write('  [../]\n')
  f_out.write('  [./min_cCr]\n')
  f_out.write('    type = NodalExtremeValue\n')
  f_out.write('    variable = cCr\n')
  f_out.write('    value_type = min\n')
  f_out.write('  [../]\n')
  f_out.write('  [./gb_cCr]\n')
  f_out.write('    type = NodalVariableValue\n')
  f_out.write('    variable = cCr\n')
  f_out.write('    nodeid = 0\n')
  f_out.write('  [../]\n')
  f_out.write('  #[./gb_cCr1]\n')
  f_out.write('  #  type = NodalVariableValue\n')
  f_out.write('  #  variable = cCr\n')
  f_out.write('  #  nodeid = 1\n')
  f_out.write('  #[../]\n')
  f_out.write('  #[./gb_cCr2]\n')
  f_out.write('  #  type = NodalVariableValue\n')
  f_out.write('  #  variable = cCr\n')
  f_out.write('  #  nodeid = 2\n')
  f_out.write('  #[../]\n')
  f_out.write('  #[./gb_cCr3]\n')
  f_out.write('  #  type = NodalVariableValue\n')
  f_out.write('  #  variable = cCr\n')
  f_out.write('  #  nodeid = 3\n')
  f_out.write('  #[../]\n')
  f_out.write('  #[./gb_cCr4]\n')
  f_out.write('  #  type = NodalVariableValue\n')
  f_out.write('  #  variable = cCr\n')
  f_out.write('  #  nodeid = 4\n')
  f_out.write('  #[../]\n')
  f_out.write('  [./max_cVa]\n')
  f_out.write('    type = NodalMaxValue\n')
  f_out.write('    variable = cVa\n')
  f_out.write('  [../]\n')
  f_out.write('  [./gb_cVa]\n')
  f_out.write('    type = NodalVariableValue\n')
  f_out.write('    variable = cVa\n')
  f_out.write('    nodeid = 0\n')
  f_out.write('  [../]\n')
  f_out.write('  [./max_cSI]\n')
  f_out.write('    type = NodalMaxValue\n')
  f_out.write('    variable = cSI\n')
  f_out.write('  [../]\n')
  f_out.write('  [./gb_cSI]\n')
  f_out.write('    type = NodalVariableValue\n')
  f_out.write('    variable = cSI\n')
  f_out.write('    nodeid = 0\n')
  f_out.write('  [../]\n')
  f_out.write('[]\n')
  f_out.write('\n')
  f_out.write('[Executioner]\n')
  f_out.write('  type = Transient\n')
  f_out.write('  #scheme = \'bdf2\'\n')
  f_out.write('  [./TimeStepper]\n')
  f_out.write('    type = SolutionTimeAdaptiveDT\n')
  f_out.write('    dt = %s # Initial time step.  In this simulation it changes.\n' %(dt))
  f_out.write('  [../]\n')
  f_out.write('  #Preconditioned JFNK (default)\n')
  f_out.write('  #solve_type = \'PJFNK\'\n')
  f_out.write('  #petsc_options_iname = \'-pc_type -pc_hypre_type -ksp_gmres_restart\'\n')
  f_out.write('  #petsc_options_value = \'hypre boomeramg 101\'\n')
  
  f_out.write('  solve_type = \'PJFNK\'\n')
  f_out.write('  petsc_options_iname = \'-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap\'\n')
  f_out.write('  petsc_options_value = \'asm         31   preonly   lu      1\'\n')
  f_out.write('\n')
  f_out.write('  l_max_its = 25\n')
  f_out.write('  nl_max_its = 15\n')
  f_out.write('  nl_rel_tol = 5.0e-6\n')
  f_out.write('  nl_abs_tol = 5.0e-7\n')
  f_out.write('  start_time = 0.0\n')
  f_out.write('  end_time = %s\n' %(end_time))
  f_out.write('  num_steps = %d\n' %(nsteps))
  f_out.write('  dtmax = 1.0e7\n')
  
  if GBbc != 'ideal':
    f_out.write('   [./TimePeriods]\n')
    f_out.write('    [./p1]\n')
    f_out.write('     start = 0.0\n')
    #f_out.write('     inactive_kernels =\'cCr_CHParsed cCr_CHint cCr_DiffusionVa cCr_DiffusionSI cVa_CHParsed cVa_Diffusion cSI_CHParsed cSI_Diffusion UISource VaSource Rate_FeCrSI Rate_FeCrVa cui_void1 cui_void2 cui_void3\'\n')
    f_out.write('     inactive_kernels =\'UISource VaSource Rate_FeCrSI Rate_FeCrVa cui_void1 cui_void2 cui_void3\'\n')
    f_out.write('    [../]\n')
    f_out.write('    [./p2]\n')
    f_out.write('      start = 1.0e-3\n')
    f_out.write('      #inactive_kernels =\'PolycrystalKernel\'\n')
    f_out.write('    [../]\n')
    f_out.write('  [../]\n')
    f_out.write('\n')
  if (GBbc != 'ideal') and False:
    f_out.write('  [./Adaptivity]\n')
    f_out.write('    print_changed_info = true\n')
    f_out.write('    initial_adaptivity = 1 #2\n')
    f_out.write('    refine_fraction = 0.3  #0.5 # Bigger = finer\n')
    f_out.write('    coarsen_fraction = 0.3 #0.2\n')
    f_out.write('    max_h_level = 20       #2\n')
    f_out.write('    weight_names = \'cCr cVa cSI gr0 gr1\'\n')
    f_out.write('    weight_values = \'1 1 1 1 1\'\n')
    f_out.write('  [../]\n')
  
  f_out.write('[]\n')
  f_out.write('\n')
  f_out.write('[Outputs]\n')
  f_out.write('  exodus = true\n')
  f_out.write('  [./csvout]\n')
  f_out.write('    type = CSV\n')
  f_out.write('    delimiter = \' \'\n')
  f_out.write('  [../]\n')
  f_out.write('  print_perf_log = true\n')
  f_out.write('[]\n')
  f_out.write('\n')

  f_out.close()
  return

# Actual Script Starts Here
if __name__ == "__main__":

    #### INPUT DATA #####
    
    casename   = 'FeCrVaSIpoly_parsed'
    
    # temp, conc, GBbc, xdim, nx, kappa_ccr
    temp_conc_list  = [ \
    ('700.0','0.03','ideal', 25.0, 125, '1.0')
    ]
    #('500.0','0.03','ideal', 25.0, 125, '1.0'),
    #('500.0','0.03','ideal', 13.9, 70, '1.0'),
    #('500.0','0.03','ideal', 9.43, 47, '1.0'),
    #('500.0','0.03','ideal', 7.46, 37, '1.0'),
    #('700.0','0.03','ideal', 25.0, 125, '1.0'),
    #('700.0','0.03','ideal', 13.9, 70, '1.0'),
    #('700.0','0.03','ideal', 9.43, 47, '1.0'),
    #('700.0','0.03','ideal', 7.46, 37, '1.0'),
    #('900.0','0.03','ideal', 25.0, 125, '1.0'),
    #('900.0','0.03','ideal', 13.9, 70, '1.0'),
    #('900.0','0.03','ideal', 9.43, 47, '1.0'),
    #('900.0','0.03','ideal', 7.46, 37, '1.0')
    
    ioffset       = 0
    machine       = 'macpro'
    nprocs        = 1 #6
    dotar         = False
    
    #xdim          = 100 #25 #25 #100
    #nx            = 150 #100 #100 #1000
    #kappa_ccr     = '100.0' #'3.5'
    
    kappa_cva     = '10.0' # <- Not Used
    kappa_csI     = '10.0' # <- Not Used
    nsteps        = 50000
    dose_rate     = 0.00001
    end_time      = '10000000.0' #1 dpa
    fthermal      = False
    use_split     = False
    
    # If 'ideal', Use Dirichlet BC's? Or use sink kernels:
    use_dirichlet = True
    
    #####################
    
    nruns = len(temp_conc_list)
    print "LENGTH of temp_conc_list = ", nruns

    for i in range(0+ioffset,nruns+ioffset):
    
      temp      = temp_conc_list[i][0]
      conc      = temp_conc_list[i][1]
      GBbc      = temp_conc_list[i][2]
      xdim      = temp_conc_list[i][3]
      nx        = temp_conc_list[i][4]
      kappa_ccr = temp_conc_list[i][5]
      print "temp, conc, GBbc = ", temp, conc, GBbc
      print "xdim, nx, kappa_ccr = ", xdim, nx, kappa_ccr
      
      if GBbc == 'ideal':
        dt = '1e-9' # '1e-5'
      else:
        dt = '1e-5' # '1e-5'
       
      if conc == '0.03':
        Ev = 1.725 # For 0.05 actually
        Ei = 3.5 # For 0.05 actually
      elif conc == '0.06':
        Ev = 1.725 # For 0.05 actually
        Ei = 3.5 # For 0.05 actually
      elif conc == '0.10':
        Ev = 1.8
        Ei = 3.5
      elif conc == '0.15':
        Ev = 1.84
        Ei = 3.5

      if temp == '300.0':
        kT = 0.00008617343 * 300.0
        if conc == '0.03':
          LCrCr_V = '1.535290E+00'
          LFeFe_V = '2.301620E+01'
          LCrFe_V = '1.922660E+00'
          LCrCr_I = '3.987870E+08'
          LFeFe_I = '3.405400E+08'
          LCrFe_I = '1.489010E+08'
        if conc == '0.06':
          LCrCr_V = '3.330430E+00'
          LFeFe_V = '3.118600E+01'
          LCrFe_V = '3.957150E+00'
          LCrCr_I = '8.863740E+08'
          LFeFe_I = '2.173130E+08'
          LCrFe_I = '1.643030E+08'
        if conc == '0.10':
          LCrCr_V = '9.414830E+00'
          LFeFe_V = '8.693120E+00'
          LCrFe_V = '2.541170E+00'
          LCrCr_I = '2.745380E+08'
          LFeFe_I = '2.071020E+08'
          LCrFe_I = '-1.225190E+08'
        if conc == '0.20':
          LCrCr_V = '2.883060E+02'
          LFeFe_V = '7.523670E+01'
          LCrFe_V = '2.883360E+01'
          LCrCr_I = '7.141370E+08'
          LFeFe_I = '4.776390E+08'
          LCrFe_I = '2.087170E+08'
      elif temp == '400.0':
        kT = 0.00008617343 * 400.0
        if conc == '0.03':
          LCrCr_V = '1.090300E+03'
          LFeFe_V = '9.348660E+03'
          LCrFe_V = '2.314570E+02'
          LCrCr_I = '3.160290E+09'
          LFeFe_I = '7.110450E+09'
          LCrFe_I = '2.332600E+08'
        if conc == '0.06':
          LCrCr_V = '1.657390E+03'
          LFeFe_V = '1.296680E+04'
          LCrFe_V = '1.449380E+03'
          LCrCr_I = '6.943490E+09'
          LFeFe_I = '3.275750E+09'
          LCrFe_I = '7.081080E+08'
        if conc == '0.10':
          LCrCr_V = '3.926240E+03'
          LFeFe_V = '1.370500E+04'
          LCrFe_V = '3.343540E+03'
          LCrCr_I = '7.377170E+09'
          LFeFe_I = '6.282770E+09'
          LCrFe_I = '3.766900E+09'
        if conc == '0.20':
          LCrCr_V = '3.131570E+04'
          LFeFe_V = '1.246480E+04'
          LCrFe_V = '-2.343070E+03'
          LCrCr_I = '4.627670E+09'
          LFeFe_I = '3.810050E+09'
          LCrFe_I = '1.985940E+09'
      elif temp == '500.0':
        kT = 0.00008617343 * 500.0
        if conc == '0.03':
          #LCrCr_V = '5.649730E+04'
          #LFeFe_V = '7.064840E+05'
          #LCrFe_V = '6.353990E+04'
          #LCrCr_I = '1.211170E+10'
          #LFeFe_I = '5.033360E+10'
          #LCrFe_I = '5.684770E+09'
          LCrCr_V = '6.311210E+04'
          LFeFe_V = '7.713660E+05'
          LCrFe_V = '1.081920E+05'
          LCrCr_I = '1.070340E+10'
          LFeFe_I = '3.219260E+10'
          LCrFe_I = '6.997280E+08'
        if conc == '0.06':
          LCrCr_V = '1.873730E+05'
          LFeFe_V = '1.025460E+06'
          LCrFe_V = '2.925230E+05'
          LCrCr_I = '2.310860E+10'
          LFeFe_I = '2.822520E+10'
          LCrFe_I = '1.388740E+09' # LCrFe_I = '-1.388740E+09'
        if conc == '0.10':
          LCrCr_V = '3.044020E+05'
          LFeFe_V = '8.132760E+05'
          LCrFe_V = '2.126840E+05'
          LCrCr_I = '4.832190E+10'
          LFeFe_I = '1.915890E+10'
          LCrFe_I = '1.731290E+10'
        if conc == '0.20':
          LCrCr_V = '3.583810E+05'
          LFeFe_V = '7.831660E+05'
          LCrFe_V = '2.587210E+05'
          LCrCr_I = '1.582010E+10'
          LFeFe_I = '2.438820E+10'
          LCrFe_I = '4.558480E+09'
      elif temp == '600.0':
        kT = 0.00008617343 * 600.0
        if conc == '0.03':
          LCrCr_V   = '7.356300E+05'
          LFeFe_V   = '7.740370E+06'
          LCrFe_V   = '1.634220E+06'
          LCrCr_I   = '3.269140E+10'
          LFeFe_I   = '1.951450E+11'
          LCrFe_I   = '3.081660E+10'
        elif conc == '0.06':
          LCrCr_V   = '1.840130E+06'
          LFeFe_V   = '9.585330E+06'
          LCrFe_V   = '1.871210E+06'
          LCrCr_I   = '5.859280E+10'
          LFeFe_I   = '1.241310E+11'
          LCrFe_I   = '2.249500E+10'
        elif conc == '0.10':
          LCrCr_V   = '4.546210E+06'
          LFeFe_V   = '1.013370E+07'
          LCrFe_V   = '2.220710E+06'
          LCrCr_I   = '1.069440E+11'
          LFeFe_I   = '1.126500E+11'
          LCrFe_I   = '8.008360E+09'
        else:
          break
      elif temp == '700.0':
        kT = 0.00008617343 * 700.0
        if conc == '0.03':
          LCrCr_V   = '4.050070E+06'
          LFeFe_V   = '7.254450E+07'
          LCrFe_V   = '4.971950E+06'
          LCrCr_I   = '4.444380E+10'
          LFeFe_I   = '4.872920E+11'
          LCrFe_I   = '2.280260E+10'
          #LCrCr_V   = '4.850020E+06'
          #LFeFe_V   = '6.168720E+07'
          #LCrFe_V   = '2.412590E+06'
          #LCrCr_I   = '4.659950E+10'
          #LFeFe_I   = '3.713600E+11'
          #LCrFe_I   = '2.064920E+10'
        elif conc == '0.06':
          LCrCr_V   = '7.961740E+06'
          LFeFe_V   = '8.508110E+07'
          LCrFe_V   = '3.628460E+06'
          LCrCr_I   = '1.729320E+11'
          LFeFe_I   = '1.694180E+11'
          LCrFe_I   = '2.608850E+10'
        elif conc == '0.10':
          LCrCr_V   = '2.046060E+07'
          LFeFe_V   = '8.943680E+07'
          LCrFe_V   = '9.753620E+06'
          LCrCr_I   = '1.541050E+11'
          LFeFe_I   = '2.576220E+11'
          LCrFe_I   = '7.116230E+10'
        else:
          break
      elif temp == '800.0':
        kT = 0.00008617343 * 800.0
        if conc == '0.03':
          LCrCr_V   = '2.286060E+07'
          LFeFe_V   = '4.042900E+08'
          LCrFe_V   = '3.686330E+07'
          LCrCr_I   = '8.894640E+10'
          LFeFe_I   = '1.163120E+12'
          LCrFe_I   = '1.240780E+11'
        elif conc == '0.06':
          LCrCr_V   = '4.674800E+07'
          LFeFe_V   = '2.699660E+08'
          LCrFe_V   = '3.212330E+07'
          LCrCr_I   = '9.904040E+10'
          LFeFe_I   = '8.759870E+11'
          LCrFe_I   = '1.363580E+11'
        elif conc == '0.10':
          LCrCr_V   = '8.956330E+07'
          LFeFe_V   = '2.674340E+08'
          LCrFe_V   = '5.677590E+07'
          LCrCr_I   = '2.788320E+11'
          LFeFe_I   = '4.537150E+11'
          LCrFe_I   = '1.537600E+11'
        else:
          break
      elif temp == '900.0':
        kT = 0.00008617343 * 900.0
        if conc == '0.03':
          #LCrCr_V   = '5.230000E+07'
          #LFeFe_V   = '8.252350E+08'
          #LCrFe_V   = '5.096550E+07'
          #LCrCr_I   = '8.693860E+10'
          #LFeFe_I   = '1.957560E+12'
          #LCrFe_I   = '4.537500E+10'
          LCrCr_V   = '6.232345E+07'
          LFeFe_V   = '6.573056E+08'
          LCrFe_V   = '2.989220E+07'
          LCrCr_I   = '8.244517E+10'
          LFeFe_I   = '1.571003E+12'
          LCrFe_I   = '1.157900E+11'
        elif conc == '0.06':
          LCrCr_V   = '8.868540E+07'
          LFeFe_V   = '7.141160E+08'
          LCrFe_V   = '-1.065620E+07'
          LCrCr_I   = '2.077440E+11'
          LFeFe_I   = '1.349830E+12'
          LCrFe_I   = '3.060640E+11'
        elif conc == '0.10':
          LCrCr_V   = '2.430980E+08'
          LFeFe_V   = '5.330930E+08'
          LCrFe_V   = '1.257910E+08'
          LCrCr_I   = '2.562710E+11'
          LFeFe_I   = '1.004860E+12'
          LCrFe_I   = '1.228320E+11'
        else:
          break
      elif temp == '1000.0':
        kT = 0.00008617343 * 1000.0
        if conc == '0.03':
          LCrCr_V   = '1.929790E+08'
          LFeFe_V   = '1.957280E+09'
          LCrFe_V   = '3.236600E+08'
          LCrCr_I   = '1.797900E+11'
          LFeFe_I   = '2.534230E+12'
          LCrFe_I   = '2.703370E+11'
        elif conc == '0.06':
          LCrCr_V   = '2.468930E+08'
          LFeFe_V   = '1.868650E+09'
          LCrFe_V   = '2.575700E+08'
          LCrCr_I   = '2.914450E+11'
          LFeFe_I   = '2.214400E+12'
          LCrFe_I   = '3.484770E+11'
        elif conc == '0.10':
          LCrCr_V   = '4.305210E+08'
          LFeFe_V   = '2.820380E+09'
          LCrFe_V   = '4.223230E+08'
          LCrCr_I   = '3.099120E+11'
          LFeFe_I   = '1.991820E+12'
          LCrFe_I   = '4.363160E+11'
        else:
          break
      else:
        break

      LCrV = -(float(LCrCr_V) + float(LCrFe_V))
      LFeV = -(float(LFeFe_V) + float(LCrFe_V))
      LCrI =  (float(LCrCr_I) + float(LCrFe_I))
      LFeI =  (float(LFeFe_I) + float(LCrFe_I))
      TestRatio = (LCrV/LFeV) - (LCrI/LFeI)
      print "TestRatio = ",TestRatio

      Xv = numpy.exp(-Ev/kT)*1.0
      Xi = numpy.exp(-Ei/kT)*1.0
      if not fthermal:
        Xv = 0.0
        Xi = 0.0
      elif Xi < 0.0: #1.0e-26:
        Xi = 0.0 #1.0e-26
 
      writeinput(casename,temp,conc,xdim,nx,LCrCr_V,LFeFe_V,LCrFe_V,LCrCr_I,LFeFe_I,LCrFe_I,nsteps,dt,GBbc,kappa_ccr,kappa_cva,kappa_csI,Xv,Xi,use_split,dose_rate,end_time,use_dirichlet)
      
      if True:
      
        if machine == 'macpro':
          os.system('mpirun -np '+str(nprocs)+' ./fecrseg/fecrseg-opt -i '+casename+'_'+GBbc+'.i >run.out')
          #os.system('peacock -i '+casename+'_'+GBbc+'.i')

        # Cleanup 1
        os.system('mkdir done.'+str(i))
        os.system('mv run.out '+casename+'* done.'+str(i))

        # Cleanup 2
        if dotar:
          os.system('tar -cvzf done.'+str(i)+'.tgz done.'+str(i))
          os.system('rm -r -f done.'+str(i))

