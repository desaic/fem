#include "StepperNewtonDyn.hpp"
#include "ElementMesh.hpp"
#include "MatrixX.hpp"
#include "ArrayUtil.hpp"
#include "femError.hpp"

#include "LinSolve.hpp"
#include "SparseLin.hpp"
#include "Eigen/Sparse"

#include <fstream>
#include <assert.h>

StepperNewtonDyn::StepperNewtonDyn():dx_tol(1e-5f),h(0.005f)
,frameCnt(0)
{
  m_Init = false;
}

void StepperNewtonDyn::init(ElementMesh * _m)
{
  m = _m;
  m->computeMass();
}

int StepperNewtonDyn::oneStep()
{
  std::vector<Vector3f> force = m->getForce();

  for(unsigned int ii =0 ; ii<force.size(); ii++){
    force[ii] = h*h*force[ii] + h*m->mass[ii]*m->v[ii];
  }

  float E = m->getEnergy();
 
  int ndof = 3*(int)m->x.size();
  std::vector<float> bb(ndof);

  compute_dx_sparse(m, force, bb);

  for(unsigned int ii =0; ii<force.size(); ii++){
    for(int jj =0 ; jj<3; jj++){
      m->x[ii][jj] += bb[3*ii+jj];
      m->v[ii][jj] = (1.0/h) * bb[3*ii+jj];
    }
  }

  //remove forces after 10 steps
  if(frameCnt==10){
    for(unsigned int ii =0; ii<force.size(); ii++){
      m->fe[ii] = Vector3f(0,0,0);
    }
  }

  frameCnt++;
  return 0;
}

float StepperNewtonDyn::compute_dx_sparse(ElementMesh * iMesh, const std::vector<Vector3f> &iForces, std::vector<float> &bb)
{
  bool triangular = true;

  if (!m_Init)
  {
    m_I.clear();
    m_J.clear();
    iMesh->stiffnessPattern(m_I, m_J, triangular);
    sparseInit();
    m_Init = true;
  }

  int ndof = bb.size();
  assert(iMesh && 3*iMesh->x.size()==bb.size());

  std::vector<float> Kvalues;
  iMesh->getStiffnessSparse(Kvalues, triangular, true);
   
  for(unsigned int ii = 0;ii<m->x.size(); ii++){
    for(int jj = 0;jj<3;jj++){
      int row = 3*ii + jj;
      if(m->fixed[ii]){
        bb[ row ] = 0;
      }else{
        bb[ row ] = iForces[ii][jj];
      }
    }
  }

  int nrow = ndof;
  int ncol = ndof;
  std::vector<int> ia(nrow+1),ja;
  std::vector<double> val ;
  std::vector<double> rhs(nrow);
  std::vector<double> x(nrow,0);
  int indVal = 0;
  //starting index of a row
  //only upper triangle is needed
  for (int col=0; col<ncol;col++){
    rhs[col] = bb[col];
    ia[col] = indVal;

    //for (Eigen::SparseMatrix<float>::InnerIterator it(K, col); it; ++it){
    //    int row = it.row();
    for (int i=m_I[col]; i<m_I[col+1]; i++)
    {
      int row = m_J[indVal];
      if(triangular && col>row)
      {
        continue;
      }
      ja.push_back(row);
      Kvalues[indVal] *= h*h;
      int vidx = row/3;
      if(row == col){
        Kvalues[indVal] += m->mass[vidx];
      }
      val.push_back(Kvalues[indVal]);
      indVal++;
    }
  }
  ia[nrow] = indVal;

  int status = sparseSolve( &(ia[0]), &(ja[0]), &(val[0]), nrow, &(x[0]), &(rhs[0]));
      std::cout<< "sparse solve " << status << "\n";
  for(int ii = 0;ii<x.size();ii++){
    bb[ii] = x[ii];
  }

  return 0;
}
