#include "FEM2DFun.hpp"
#include "ElementMesh2D.h"
#include "RealField.hpp"

void FEM2DFun::setParam(const Eigen::VectorXd & x0)
{

}

double FEM2DFun::f()
{
  return 0;
}

Eigen::VectorXd FEM2DFun::df()
{
  return Eigen::VectorXd(1);
}

FEM2DFun::FEM2DFun() :em(0), field(0)
{
}

FEM2DFun::~FEM2DFun(){}