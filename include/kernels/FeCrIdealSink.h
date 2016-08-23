#ifndef FeCrIdealSink_H
#define FeCrIdealSink_H

#include "Kernel.h"
#include "Material.h"

//Forward Declarations
class FeCrIdealSink;

template<>
InputParameters validParams<FeCrIdealSink>();

class FeCrIdealSink : public Kernel
{
public:

  FeCrIdealSink(const InputParameters & parameters);
  
protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();

private:
  
  //std::vector<VariableValue *> _vals;
  //std::vector<VariableGradient *> _grad_vals;
  unsigned int _n;
  const VariableValue & _c;
  Real _xdim;
  Real _dose_rate;
  Real _fcut;
  std::string _diffusivitySink;
  const MaterialProperty<Real> & _diffusivity;

};
#endif //FeCrIdealSink

