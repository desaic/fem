#ifndef STEPPERNEWTON_HPP
#define STEPPERNEWTON_HPP
#include "Stepper.hpp"
class StepperNewton:public Stepper
{
public:
  StepperNewton();
  int oneStep();
  float oneStepDense();
  float oneStepSparse();
  bool dense;

  float dx_tol;
  float h;
  bool rmRigid;
};
#endif
