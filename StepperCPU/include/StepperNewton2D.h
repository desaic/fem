#ifndef STEPPERNEWTON2D_HPP
#define STEPPERNEWTON2D_HPP
#include "Stepper2D.h"
#include <vector>
#include "cfgDefs.h"

class StepperNewton2D:public Stepper2D
{
public:
  StepperNewton2D();
  int oneStep();

  void enforcePeriodicity(bool iOn);
  void removeTranslation(bool iOn);
  void removeRotation(bool iOn);

  bool dense;
  cfgScalar dx_tol;
  cfgScalar h;

private:
  int compute_dx_dense(ElementMesh2D * iMesh, const std::vector<Vector2S> &iForces, bool iRmRigid, std::vector<cfgScalar> &bb);
  int compute_dx_sparse(ElementMesh2D * iMesh, const std::vector<Vector2S> &iForces, bool iRemoveTranslation, bool iRemoveRotation, bool iPeriodic, std::vector<cfgScalar> &bb);

private:
  std::vector<int> m_I, m_J;
  bool m_Init;

  bool m_periodic;
  bool m_noTranslation;
  bool m_noRotation;
};
#endif
