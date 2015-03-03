#ifndef NEWTON_CUDA_HIER_HPP
#define NEWTON_CUDA_HIER_HPP
#include "Stepper.hpp"
#include "LinSolveCusp.hpp"
class NewtonCudaHier:public Stepper
{
public:
  NewtonCudaHier();
  void init(ElementMesh * _m);
  int oneStep();
  LinSolveCusp solver;

  float dx_tol;
  float h;
  //number of iterations for linear solver;
  int linIter;
  //level it is solving
  int level;
};

#endif