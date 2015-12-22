#ifndef STEPPERNEWTONDYN_HPP
#define STEPPERNEWTONDYN_HPP
#include "Stepper.hpp"
#include <vector>
#include <Eigen/Dense>
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
  ///@param collide. Pass in empty array if no collision
  float compute_dx_sparse(ElementMesh * iMesh, const std::vector<Eigen::Vector3f> &iForces,
                          const std::vector<bool> &collide,
                          std::vector<float> &bb);

private:
  std::vector<int> m_I, m_J;
  bool m_Init;
};
#endif
