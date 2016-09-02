#include "StepperGrad.hpp"
#include "ElementMesh.hpp"
#include "femError.hpp"
#include "ArrayUtil.hpp"

#include <vector>
#include <iostream>
using namespace Eigen;
StepperGrad::StepperGrad():h(0.01f),force_L2tol(1e-4f) {}

int StepperGrad::oneStep()
{
  std::vector<Vector3S> force = m->getForce();
  for(unsigned int ii = 0;ii<force.size();ii++){
    if(m->fixed[ii]){
      force[ii] = Vector3S(0,0,0);
    }
  }
  cfgScalar E = m->getEnergy();
  cfgScalar totalMag = 0;
  for(unsigned int ii = 0;ii<force.size();ii++){
    totalMag += force[ii].squaredNorm();  
  }
  if(totalMag<force_L2tol){
    return 0;
  }

  std::vector<Vector3S> x0 = m->x;
  while(1){
    m->x=x0;
    addmul(m->x, h, force);
    cfgScalar E1 = m->getEnergy();

    if(E1>E || fem_error){
      fem_error = 0;
      h = 0.5f* h;
      std::cout<<"h "<<h<<"\n";
    }else{
      h=1.1f*h;
      break;
    }
  }
  return 0;
}
