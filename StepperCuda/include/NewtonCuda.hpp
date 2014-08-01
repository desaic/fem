#ifndef NEWTONCUDA_HPP
#define NEWTONCUDA_HPP
#include "Stepper.hpp"

class ElementMesh;
class ConjugateGradientCuda;
class NewtonCuda:public Stepper
{
public:
  NewtonCuda();
  ///must be called to initialize or if topology of the mesh changes
  int initCuda(ElementMesh * _m);
  void step(ElementMesh * _m=0);

  ConjugateGradientCuda * solver;
  ElementMesh * m;
};
#endif