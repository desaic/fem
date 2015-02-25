#ifndef STEPPER_CUDA_HPP
#define STEPPER_CUDA_HPP
#include "Stepper.hpp"
#include "LinSolveCusp.hpp"
class NewtonCuda :public Stepper
{
public:
  NewtonCuda();
  void init(ElementMesh * _m);
  int oneStep() ;
  LinSolveCusp solver;

  float dx_tol;
  float h;
};

#endif