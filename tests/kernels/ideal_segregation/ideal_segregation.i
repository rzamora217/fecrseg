# FeCrSeg Test Input File: ideal_segregation.i 

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 37
  ny = 1
  nz = 0
  xmin = 0
  xmax = 7
  ymin = 0
  ymax = 1
  zmin = 0
  zmax = 0
  elem_type = QUAD4
[]

[Variables]
  [./cCr]
    order = THIRD
    family = HERMITE
    scaling=1.00E-00
    block = 0
  [../]
  [./cVa]
    order = THIRD #FIRST
    family = HERMITE #LAGRANGE
    scaling=1.00E-00
    block = 0
  [../]
  [./cSI]
    order = THIRD #FIRST
    family = HERMITE #LAGRANGE
    scaling=1.00E-00
    block = 0
  [../]
[]

[Materials]
  active = 'LANL free_energy'
  [./LANL]

    type = FeCrVaSIbulk

    cCr = cCr
    cVa = cVa
    cSI = cSI

    kappa_cfe = 1.0
    kappa_cva = 10.0
    kappa_ccr = 1.0
    kappa_csI = 10.0

    ICrva = 0.0 #-0.048 # eV
    ICrsI = 0.0 #-0.065 # eV
    IVava = 0.0 #-0.049 # eV

    temp = 500.0 # K
    Laav1 = 6.311210E+04
    Lbbv1 = 7.713660E+05
    Labv1 = 1.081920E+05
    Laai1 = 1.070340E+10
    Lbbi1 = 3.219260E+10
    Labi1 = 6.997280E+08

    Z = 1.0
    dHVaUI = 1.2 #eV
    block = 0
  [../]
  [./free_energy]
    type = DerivativeParsedMaterial
    block = 0
    f_name = F
    args = 'cCr cVa cSI'
    constant_names       = 'T0  y   kB'
    constant_expressions = '410 500.0 8.6173324e-5'
    derivative_order = 2
    enable_jit = true

    # G(CFe) from Mathematica
    function = '(cCr*-4.16723234898182 + (1 - cCr)*-3.85953051397515 + cCr*(1 - cCr)*(0.380577499303532 + 0.0885190880468175*(1 - 2*cCr) + -0.0549325534161349*(1 - 2*cCr)^2 + 0.203752174307515*(1 - 2*cCr)^3 + -0.164028794268281*(1 - 2*cCr)^4 + -0.0172550241855144*(1 - 2*cCr)^5)) * y/T0
    + (cCr*-4.12297413272388 + (1 - cCr)*-3.83660056975999 + cCr*(1 - cCr)*(0.385512095260579 + 0.0962430189496054*(1 - 2*cCr) + -0.041249704177587*(1 - 2*cCr)^2 + 0.194439246552959*(1 - 2*cCr)^3 + -0.195412847295217*(1 - 2*cCr)^4 + 0.00967038578662529*(1 - 2*cCr)^5))
    -(0.000263056717498972 + 4.614980081531*10^-5*cCr + -4.75235048526914*10^-5*cCr^2 + 9.35759929354588*10^-6*cCr^3)*y*log(y)
    - (3.04676203180853*10^-9 + -2.07225774483557*10^-8*cCr + 3.55582178830517*10^-8*cCr^2 + -2.70425743485173*10^-8*cCr^3)*y^2
    + (-1.51823088659839*10^-13 + 5.18553402098699*10^-12*cCr + -4.56309143596694*10^-12*cCr^2 + 1.08597105154957*10^-11*cCr^3)/2*y^3
    + (-(cCr*-4.12297413272388 + (1 - cCr)*-3.83660056975999 +
     cCr*(1 - cCr)*(0.385512095260579 +
         0.0962430189496054*(1 - 2*cCr) + -0.041249704177587*(1 - 2*cCr)^2 +
             0.194439246552959*(1 - 2*cCr)^3 + -0.195412847295217*(1 - 2*cCr)^4 +
                 0.00967038578662529*(1 - 2*cCr)^5))/T0 + (3.04676203180853*10^-9 + -2.07225774483557*10^-8*cCr +
                  3.55582178830517*10^-8*cCr^2 + -2.70425743485173*10^-8*cCr^3)*T0 + (0.000263056717498972 +
                   4.614980081531*10^-5*cCr + -4.75235048526914*10^-5*cCr^2 +
                    9.35759929354588*10^-6*cCr^3)*log(T0) + (-1.51823088659839*10^-13 +
                     5.18553402098699*10^-12*cCr + -4.56309143596694*10^-12*cCr^2 +
                      1.08597105154957*10^-11*cCr^3)/2*T0^2)*y + y*kB*(cCr*plog(cCr,0.0001) + (1 - cCr)*plog(1 - cCr,0.0001))'
  [../]
[]

[Kernels]
  active = 'cCr_CHParsed cCr_CHint cCr_DiffusionVa cCr_DiffusionSI ie_cCr cVa_CHParsed cVa_Diffusion cVa_CHint ie_cVa cSI_CHParsed cSI_Diffusion cSI_CHint ie_cSI UISource VaSource Rate_FeCrSI Rate_FeCrVa'
  [./Va_idealsink]
    type = FeCrIdealSink
    variable = cVa
    c = cSI
    xdim = 7.460000
    dose_rate = 0.000010
    fcut = 0.100000
    diffusivitySink = DrecombVa
  [../]
  [./SI_idealsink]
    type = FeCrIdealSink
    variable = cSI
    c = cVa
    xdim = 7.460000
    dose_rate = 0.000010
    fcut = 0.100000
    diffusivitySink = DrecombSI
  [../]
  [./Cr_GBsink]
    type = CHRadSinkFeCrNoPF
    xdim = 7.460000
    variable = cCr
    mob_name = M_Cr
    LogC_name = LogC
    LogTol_name = LogTol
  [../]
  [./Va_GBsink]
    type = CHRadSinkFeVaNoPF
    xdim = 7.460000
    variable = cVa
    mob_name = M_Va_sink
    LogC_name = LogC
    LogTol_name = LogTol
  [../]
  [./SI_GBsink]
    type = CHRadSinkFeSINoPF
    xdim = 7.460000
    variable = cSI
    mob_name = M_SI_sink
    LogC_name = LogC
    LogTol_name = LogTol
  [../]
  [./cCr_CHRes]
    type = SplitCHParsed
    variable = cCr
    f_name = F
    kappa_name = kappa_cCr
    w = wCr
  [../]
  [./wCr_CHRes]
    type = SplitCHWRes
    variable = wCr
    mob_name = M_Cr
  [../]
  [./time_cCr]
    type = CoupledTimeDerivative
    variable = wCr
    v = cCr
  [../]

  [./cCr_CHParsed]
    type = CahnHilliard
    variable = cCr
    f_name = F
    #args = cCr
    mob_name = M_Cr
  [../]
  [./cCr_CHint]
    type = CHInterface
    variable = cCr
    mob_name = M_Cr
    kappa_name = kappa_cCr
  [../]
  [./cCr_DiffusionVa]
    type = DiffusionAdd
    variable = cCr
    diffusivity = dCrV
    c = cVa
    dirscale = -1.0
  [../]
  [./cCr_DiffusionSI]
    type = DiffusionAdd
    variable = cCr
    diffusivity = dCrI
    c = cSI
    dirscale = 1.0
  [../]
  [./ie_cCr]
    type = TimeDerivative
    variable = cCr
  [../]

  [./cVa_CHParsed]
    type = CahnHilliardD
    variable = cVa
    f_name = F
    args = cCr
    mob_name = M_Va
  [../]
  [./cVa_CHint]
    type = CHInterface
    variable = cVa
    mob_name = M_Va_int
    kappa_name = kappa_cVa
  [../]
  [./cVa_Diffusion]
    type = DiffusionAdd2
    variable = cVa
    diffusivityCr = dCrV
    diffusivityFe = dFeV
    cCr = cCr
    cVa = cVa
    xdim = 7.460000
    fcut = 0.100000
  [../]
  [./ie_cVa]
    type = TimeDerivative
    variable = cVa
  [../]

  [./cSI_CHParsed]
    type = CahnHilliardD
    variable = cSI
    f_name = F
    args = cCr
    mob_name = M_SI
  [../]
  [./cSI_CHint]
    type = CHInterface
    variable = cSI
    mob_name = M_SI_int
    kappa_name = kappa_cCr
  [../]
  [./cSI_Diffusion]
    type = DiffusionAdd2
    variable = cSI
    diffusivityCr = dCrI
    diffusivityFe = dFeI
    cCr = cCr
    cVa = cVa
    xdim = 7.460000
    fcut = 0.100000
  [../]
  [./ie_cSI]
    type = TimeDerivative
    variable = cSI
  [../]

  [./Rate_FeCrSI]
    type = FeCrDRecombRate
    variable = cSI
    diffusivityCrI = dCrI
    diffusivityFeI = dFeI
    diffusivityCrV = dCrV
    diffusivityFeV = dFeV
    cCr = cCr
    c = cVa
    cVa = cVa
    cSI = cSI
    LogC_name = LogC
    LogTol_name = LogTol
  [../]
  [./Rate_FeCrVa]
   type = FeCrDRecombRate
    variable = cVa
    diffusivityCrI = dCrI
    diffusivityFeI = dFeI
    diffusivityCrV = dCrV
    diffusivityFeV = dFeV
    cCr = cCr
    c = cSI
    cVa = cVa
    cSI = cSI
    LogC_name = LogC
    LogTol_name = LogTol
  [../]

  [./UISource]
    type = UISource
    variable = cSI
    R = 0.000010
    Factor = 1.0
  [../]
  [./VaSource]
    type = VaSource
    variable = cVa
    R = 0.000010
    Factor = 1.0
  [../]
[]

[BCs]
  active = 'cVa_gb cSI_gb Periodic'
  [./cVa_gb]
    type = DirichletBC
    variable = cVa
    boundary = 'left right'
    value = 0.000000e+00
  [../]
  [./cVa_gbn]
    type = NeumannBC
    variable = cVa
    boundary = 'left right'
    value = 0.0
  [../]
  [./cSI_gb]
    type = DirichletBC
    variable = cSI
    boundary = 'left right'
    value = 0.000000e+00
  [../]
  [./cSI_gbn]
    type = NeumannBC
    variable = cSI
    boundary = 'left right'
    value = 0.0
  [../]
  [./cCr_gbl]
    type = NeumannBC
    variable = cCr
    boundary = 'left'
    value = 0.0
  [../]
  [./cCr_gbr]
    type = NeumannBC
    variable = cCr
    boundary = 'right'
    value = 0.0
  [../]
  [./Periodic]
    [./cCr]
      variable = cCr
      auto_direction = 'x y'
    [../]
    [./cVa]
      variable = cVa
      auto_direction = 'x y'
    [../]
    [./cSI]
      variable = cSI
      auto_direction = 'x y'
    [../]
  [../]
[]


[ICs]
  active = 'cCr cVa cSI'
  [./cCr]
    type = ConstantIC
    variable = cCr
    value = 0.03
  [../]
  [./cVa]
    type = ConstantIC
    variable = cVa
    value = 0.000000e+00
  [../]
  [./cSI]
    type = ConstantIC
    variable = cSI
    value = 0.000000e+00
  [../]
  [./cCr_Random]
    type = RandomIC
    variable = cCr
    min = 0.05
    max = 0.15
    block = 0
  [../]
  [./cCr_Point]
    type = BoundingBoxIC
    variable = cCr
    x1 = -0.001 # 4.999 #-0.001 # x coordinate of the lower left-hand corner of the box
    y1 = -0.001 #-0.001 #-0.001 # y coordinate of the lower left-hand corner of the box
    x2 =  0.001 # 5.001 # 0.001 # x coordinate of the upper right-hand corner of the box
    y2 =  1.001 # 1.001 # 1.001 # y coordinate of the upper right-hand corner of the box
    inside = 0.1
    outside = 0.03
    block = 0
  [../]
[]


[Executioner]
  type = Transient
  #scheme = 'bdf2'
  [./TimeStepper]
    type = SolutionTimeAdaptiveDT
    dt = 1e-9 # Initial time step.  In this simulation it changes.
  [../]
  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -ksp_grmres_restart -sub_ksp_type -sub_pc_type -pc_asm_overlap'
  petsc_options_value = 'asm         31   preonly   lu      1'

  l_max_its = 25
  nl_max_its = 15
  nl_rel_tol = 5.0e-6
  nl_abs_tol = 5.0e-7
  start_time = 0.0
  end_time = 10000.0
  num_steps = 100
  dtmax = 1.0e7
[]

[Outputs]
  exodus = true
[]

