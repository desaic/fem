#include "StepperNewton.hpp"
#include "ElementMesh.hpp"
#include "MatrixXd.hpp"
#include <complex>
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include "lapacke.h"

StepperNewton::StepperNewton():dense(true)
{}

float StepperNewton::oneStepSparse(ElementMesh * m)
{
  return -1;
}

float StepperNewton::oneStepDense(ElementMesh * m)
{
  std::vector<Vector3f> force = m->getForce();
  float E0 = m->getEnergy();
  MatrixXd K = m->getStiffness();
  return E0;
}

void StepperNewton::step(ElementMesh * m)
{
  for(int ii = 0;ii<nSteps;ii++){
    float E =0;
    if(dense){
      E = oneStepDense(m);
    }else{
      E = oneStepSparse(m);
    }
  }
}