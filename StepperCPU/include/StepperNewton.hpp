#ifndef STEPPERNEWTON_HPP
#define STEPPERNEWTON_HPP
#include "Stepper.hpp"
class StepperNewton:public Stepper
{
public:
  StepperNewton();
  void step(ElementMesh * m);
  float oneStepDense(ElementMesh * m);
  float oneStepSparse(ElementMesh * m);
  bool dense;

  double * bb ;
  float force_L2tol;
  float h;
};
#endif