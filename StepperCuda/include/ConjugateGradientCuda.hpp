#ifndef CONJUGATEGRADIENTCUDA_HPP
#define CONJUGATEGRADIENTCUDA_HPP
class ConjugateGradientCuda
{
public:
  ConjugateGradientCuda();
  int initCuda();
  int solve(int N, int nz , int * I, int * J, float * val,
                          float * x, float * rhs);
  float tol;
  int devID;
};
#endif