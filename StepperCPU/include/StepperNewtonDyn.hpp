#ifndef STEPPERNEWTONDYN_HPP
#define STEPPERNEWTONDYN_HPP
#include "Stepper.hpp"
#include <vector>
#include "Vector3f.h"

class StepperNewtonDyn:public Stepper
{
public:
  StepperNewtonDyn();
  virtual void init(ElementMesh * _m);
  int oneStep();

  float dx_tol;
  float h;

  int frameCnt;
private:
  float compute_dx_sparse(ElementMesh * iMesh, const std::vector<Vector3f> &iForces, std::vector<float> &bb);

private:
  std::vector<int> m_I, m_J;
  bool m_Init;
};
#endif
