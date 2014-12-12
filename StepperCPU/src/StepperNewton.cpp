#include "StepperNewton.hpp"
#include "ElementMesh.hpp"
#include "MatrixXd.hpp"
#include "ArrayUtil.hpp"
#include "femError.hpp"

#include "LinSolve.hpp"

StepperNewton::StepperNewton():dense(true),dx_tol(1e-5f),h(1.0f)
,rmRigid(false)
{}

float StepperNewton::oneStepSparse(ElementMesh * m)
{
  return -1;
}

///@brief add rows to K and b to constrain 6 degrees of freedom.
///@param K size is #DOF + 6
void fixRigid(MatrixXd & K, double * b,
        ElementMesh * mesh)
{
  int row = K.mm-6;
  for(int ii =0; ii<6; ii++){
    b[row + ii] = 0;
  }

  Vector3f center(0,0,0);
  for(int ii = 0;ii<mesh->x.size();ii++){
    center += mesh->x[ii];
  }
  center /= (float)mesh->x.size();
  float cscale = 1000;
  for(int ii = 0;ii<mesh->x.size();ii++){
    int col = 3*ii;
    //relative position to center
    Vector3f rel = mesh->x[ii] - center;
    //cross product matrix
    Matrix3f c (0, -rel[2],  rel[1],
             rel[2],       0, -rel[0],
            -rel[1],  rel[0],       0);

    c = cscale * c;
    for(int jj=0; jj<3; jj++){
      for(int kk = 0; kk<3; kk++){
        K(row+jj, col+kk) = c(jj,kk);
        K(col+kk, row+jj) = c(jj,kk);
      }
    }

    for(int jj=0; jj<3; jj++){
      K(row+3+jj, col+jj) = 1;
      K(col+jj, row+3+jj) = 1;
    }
  }
}

float StepperNewton::oneStepDense(ElementMesh * m)
{
  std::vector<Vector3f> force = m->getForce();
  float E = m->getEnergy();

  MatrixXd K = m->getStiffness();
  int ndof = 3*(int)m->x.size();
  std::vector<double> bb (ndof);

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

  if(rmRigid){
    K.resize(ndof + 6, ndof+6);
    bb.resize(ndof+6);
    fixRigid(K,&(bb[0]),m);
  }

  linSolve(K,&(bb[0]));

  double totalMag = 0;
  for(int ii = 0;ii<ndof;ii++){
    totalMag += std::abs(bb[ii]);
  }

  for(unsigned int ii = 0;ii<m->x.size(); ii++){
    for(int jj = 0;jj<3;jj++){
      force[ii][jj] = (float)bb[3*ii+jj];
    }
  }
  //line search
  std::vector<Vector3f> x0 = m->x;
  float E1=E;
  while(1){
    if(totalMag * h<dx_tol){
      break;
    }
    m->x=x0;
    addmul(m->x, h, force);
    E1 = m->getEnergy();
    if(E1>E || fem_error){
      fem_error = 0;
      h = 0.5f* h;
      if(h<1e-12){
        break;
      }
    }else{
      h=1.1f*h;
      break;
    }
  }

  return E1;
}

void StepperNewton::step(ElementMesh * m)
{

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
}
