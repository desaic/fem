#include "NewtonCuda.hpp"
#include "ConjugateGradientCuda.hpp"
NewtonCuda::NewtonCuda():solver(0),m(0)
{
  solver = new ConjugateGradientCuda();
}

int NewtonCuda::initCuda(ElementMesh * _m)
{
  m=_m;
  int N, nz;
  int * I, * J;
  solver->initCuda(N,nz,I,J);
  return 0;
}