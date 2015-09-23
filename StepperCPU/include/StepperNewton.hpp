#ifndef STEPPERNEWTON_HPP
#define STEPPERNEWTON_HPP
#include "Stepper.hpp"
#include <vector>
#include "Vector3f.h"

class StepperNewton:public Stepper
{
public:
  StepperNewton();
  int oneStep();
  bool dense;

  float dx_tol;
  float h;
  bool rmRigid;

private:
  int compute_dx_dense(ElementMesh * iMesh, const std::vector<Vector3f> &iForces, bool iRmRigid, std::vector<float> &bb);
  int compute_dx_sparse(ElementMesh * iMesh, const std::vector<Vector3f> &iForces, bool iRmRigid, std::vector<float> &bb);

private:
  std::vector<int> m_I, m_J;
  bool m_Init;
};
#endif
