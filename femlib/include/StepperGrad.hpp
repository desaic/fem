#ifndef STEPPERGRAD_HPP
#define STEPPERGRAD_HPP
#include "Stepper.hpp"
#include "cfgDefs.h"
///@brief very slow example time stepper
class StepperGrad:public Stepper
{
public:
  StepperGrad();
  int oneStep();
  cfgScalar h;
  ///@brief tolerance for sum of force length squared.
  cfgScalar force_L2tol;
};
#endif