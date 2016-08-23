/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "CahnHilliardD.h"

template<>
InputParameters validParams<CahnHilliardD>()
{
  InputParameters params = CahnHilliardBaseD<Real>::validParams();
  params.addClassDescription("Cahn-Hilliard Kernel that uses a DerivativeMaterial Free Energy and a scalar (isotropic) mobility");
  return params;
}

CahnHilliardD::CahnHilliardD(const InputParameters & parameters) :
    CahnHilliardBaseD<Real>(parameters)
{
}

