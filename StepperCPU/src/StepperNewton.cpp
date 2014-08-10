#include "StepperNewton.hpp"
#include "ElementMesh.hpp"
#include "MatrixXd.hpp"
#include "ArrayUtil.hpp"
#include "femError.hpp"

#include "LinSolve.hpp"

StepperNewton::StepperNewton():dense(true),bb(0),force_L2tol(1e-1f),h(1.0f)
{}

float StepperNewton::oneStepSparse(ElementMesh * m)
{
  return -1;
}

float StepperNewton::oneStepDense(ElementMesh * m)
{
  std::vector<Vector3f> force = m->getForce();
  float E = m->getEnergy();
  float totalMag = 0;
  for(unsigned int ii = 0;ii<force.size();ii++){
    totalMag += force[ii].absSquared();  
  }
  if(totalMag<force_L2tol){
    return E;
  }

  MatrixXd K = m->getStiffness();
  
  int ndof = 3*(int)m->x.size();
  for(unsigned int ii = 0;ii<m->x.size(); ii++){
    for(int jj = 0;jj<3;jj++){
      int row = 3*ii + jj;
      K(row,row) += 100;
      if(m->fixed[ii]){
        bb[ row ] = 0;
        for(int kk = 0;kk<ndof;kk++){
          K(row,kk) = 0;
          K(row,row) = 1;
        }
      }else{
        bb[ row ] = force[ii][jj];
      }
    }
  }
  linSolve(K,bb);

  for(unsigned int ii = 0;ii<m->x.size(); ii++){
    for(int jj = 0;jj<3;jj++){
      force[ii][jj] = (float)bb[3*ii+jj];
    }
  }

  //line search
  std::vector<Vector3f> x0 = m->x;
  float E1;
  while(1){
    m->x=x0;
    addmul(m->x, h, force);
    E1 = m->getEnergy();

    if(E1>E || fem_error){
      fem_error = 0;
      h = 0.5f* h;
      std::cout<<"h "<<h<<"\n";
    }else{
      h=1.1f*h;
      break;
    }
  }
  return E1;
}

void StepperNewton::step(ElementMesh * m)
{
  int ndof = 3*(int)m->x.size();
  bb = new double[ndof];
  float E0 = m->getEnergy();
  for(int ii = 0;ii<nSteps;ii++){
    float E =0;
    if(dense){
      E = oneStepDense(m);
    }else{
      E = oneStepSparse(m);
    }
    if(std::abs(E0-E)<1e-3f){
      break;
    }
    std::cout<<E<<"\n";
    E0=E;
  }
  
  delete bb;
}