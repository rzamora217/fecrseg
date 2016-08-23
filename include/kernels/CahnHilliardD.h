/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CAHNHILLIARDD_H
#define CAHNHILLIARDD_H

#include "CahnHilliardBaseD.h"

/**
 * SplitCHWRes creates the residual of the Cahn-Hilliard
 * equation with a scalar (isotropic) mobility.
 */
class CahnHilliardD : public CahnHilliardBaseD<Real>
{
public:
  CahnHilliardD(const InputParameters & parameters);
};

template<>
InputParameters validParams<CahnHilliardD>();

#endif // CAHNHILLIARDD_H
