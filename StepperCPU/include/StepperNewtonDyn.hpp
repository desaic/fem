#ifndef STEPPERNEWTONDYN_HPP
#define STEPPERNEWTONDYN_HPP
#include "Stepper.hpp"
#include <vector>
#include <Eigen/Dense>
#include "cfgDefs.h"
class StepperNewtonDyn:public Stepper
{
public:
  StepperNewtonDyn();
  virtual void init(ElementMesh * _m);
  int oneStep();

  cfgScalar dx_tol;
  cfgScalar h;

  int frameCnt;
private:
  ///@param collide. Pass in empty array if no collision
  cfgScalar compute_dx_sparse(ElementMesh * iMesh, const std::vector<Vector3S> &iForces,
                          const std::vector<bool> &collide, 
                          std::vector<cfgScalar> &bb);

private:
  std::vector<int> m_I, m_J;
  bool m_Init;
};
#endif
