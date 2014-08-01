#include "StepperGrad.hpp"
#include "ElementMesh.hpp"
#include "femError.hpp"
#include "ArrayUtil.hpp"

#include <vector>
#include <iostream>
StepperGrad::StepperGrad():h(0.01f),force_L2tol(1e-4f) {}

///@return energy value.
float StepperGrad::oneStep(ElementMesh * m)
{
  std::vector<Vector3f> force = m->getForce();
  for(unsigned int ii = 0;ii<force.size();ii++){
    if(m->fixed[ii]){
      force[ii] = Vector3f(0,0,0);
    }
  }
  float E = m->getEnergy();
  float totalMag = 0;
  for(unsigned int ii = 0;ii<force.size();ii++){
    totalMag += force[ii].absSquared();  
  }
  if(totalMag<force_L2tol){
    return E;
  }

  std::vector<Vector3f> x0 = m->x;
  while(1){
    m->x=x0;
    addmul(m->x, h, force);
    float E1 = m->getEnergy();

    if(E1>E || fem_error){
      fem_error = 0;
      h = 0.5f* h;
      std::cout<<"h "<<h<<"\n";
    }else{
      h=1.1f*h;
      break;
    }
  }
  return E;
}

void StepperGrad::step(ElementMesh * m)
{
  for(int ii = 0;ii<nSteps;ii++){
    float E = oneStep(m);
    std::cout<<"E "<<E<<"\n";
  }
}