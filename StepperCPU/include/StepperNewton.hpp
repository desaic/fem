#ifndef STEPPERNEWTON_HPP
#define STEPPERNEWTON_HPP
#include "Stepper.hpp"
#include <vector>
#include <Eigen/Dense>
#include "cfgDefs.h"
class StepperNewton:public Stepper
{
public:
  StepperNewton();
  int oneStep();

  void enforcePeriodicity(bool iOn);
  void removeTranslation(bool iOn);
  void removeRotation(bool iOn);

  bool dense;
  cfgScalar dx_tol;
  cfgScalar h;

private:
  int compute_dx_dense(ElementMesh * iMesh, const std::vector<Vector3S> &iForces, bool iRmRigid, std::vector<cfgScalar> &bb);
  int compute_dx_sparse(ElementMesh * iMesh, const std::vector<Vector3S> &iForces, bool iRmRigid, std::vector<cfgScalar> &bb);
  int compute_dx_sparse(ElementMesh * iMesh, const std::vector<Vector3S> &iForces, bool iRemoveTranslation, bool iRemoveRotation, bool iPeriodic, std::vector<cfgScalar> &bb);

private:
  std::vector<int> m_I, m_J;
  bool m_Init;

  bool m_periodic;
  bool m_noTranslation;
  bool m_noRotation;
};
#endif
