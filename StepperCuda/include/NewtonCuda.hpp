#ifndef STEPPER_CUDA_HPP
#define STEPPER_CUDA_HPP
#include "Stepper.hpp"
class NewtonCuda :public Stepper
{
public:
  NewtonCuda();
  virtual void init(ElementMesh * _m);
  virtual int oneStep() = 0;
};

#endif