#include "StepperNewton2D.h"
#include "ElementMesh2D.h"
#include "ArrayUtil.hpp"
#include "femError.hpp"
#include "MeshUtilities.h"

#include "SparseLin.hpp"
#include "Eigen/Sparse"

#include <fstream>
#include <assert.h>
#include <iostream>

StepperNewton2D::StepperNewton2D():dense(true),dx_tol(1e-5f),h(1.0f)
{
  m_Init = false;
  dense = false;

  m_periodic = false;
  m_noTranslation = false;
  m_noRotation = false;
}

void StepperNewton2D::enforcePeriodicity(bool iOn)
{
  m_periodic = iOn;
}

void StepperNewton2D::removeTranslation(bool iOn)
{
  m_noTranslation = iOn;
}

void StepperNewton2D::removeRotation(bool iOn)
{
  m_noRotation = iOn;
}

int StepperNewton2D::oneStep()
{
  int status = 0;

  fem_error = FEM_OK;

  std::vector<Vector2S> force = m->getForce();
  cfgScalar E = m->getEnergy();
 
  int ndof = 2*(int)m->x.size();
  int nconstraint = 0;
  if (m_noTranslation)
  {
    nconstraint +=2;
  }
  if (m_noRotation)
  {
    nconstraint +=1;
  }
  if (m_periodic)
  {
    std::vector<int> sideVertexIndices[4];
    int iside;
    for (iside=0; iside<4; iside++)
    {
        meshUtil::getSideVertices(iside, (const ElementRegGrid2D*)m, sideVertexIndices[iside]);
    }
    nconstraint += 2*((int)sideVertexIndices[1].size()-1 + (int)sideVertexIndices[3].size()-2);
  }
  std::vector<cfgScalar> bb(ndof+nconstraint);

  if (dense)
  {
    assert(0); // NOT SUPPORTED FOR NOW
    //status = compute_dx_dense(m, force, rmRigid, bb);
  }
  else
  {
    status = compute_dx_sparse(m, force, m_noTranslation, m_noRotation, m_periodic, bb);
  }
  if (status>0)
    return status;

  double maxVal = -FLT_MAX;
  double totalMag = 0;
  for(int ii = 0;ii<ndof;ii++){
    totalMag += std::abs(bb[ii]);
    /*if (abs(bb[ii])>maxVal)
    {
      maxVal = abs(bb[ii]);
    }*/ 
  }
  //std::cout << "norm forces = " << totalMag << ", max val = " << maxVal << std::endl;

  for(unsigned int ii = 0;ii<m->x.size(); ii++){
    for(int jj = 0;jj<2;jj++){
      force[ii][jj] = (cfgScalar)bb[2*ii+jj];
    }
  }
  //line search
  h = 1;
  std::vector<Vector2S> x0 = m->x;
  cfgScalar E1=E;
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
     // h=1.1f*h;
      break;
    }
  }
  //std::cout << "h = " << h << std::endl;
  if (std::abs(E - E1) < 1e-3f){
    return -1;
  }
  //std::cout << E << "\n";
  return 0;
}

///@brief add rows to K and b to constrain 3 degrees of freedom.
///@param K size is #DOF + 3
void fixRigid(MatrixXS & K, cfgScalar * b,  ElementMesh2D * mesh)
{
  int row = K.rows()-3;
  for(int ii =0; ii<3; ii++){
    b[row + ii] = 0;
  }

  Vector2S center(0,0);
  for(int ii = 0;ii<mesh->x.size();ii++){
    center += mesh->x[ii];
  }
  center /= (cfgScalar)mesh->x.size();
  cfgScalar cscale = 1000;
  for(int ii = 0;ii<mesh->x.size();ii++){
    int col = 2*ii;
    //relative position to center
    Vector2S rel = mesh->x[ii] - center;
    //cross product matrix
  /*  Vector2S c (0, -rel[2],  rel[1],
             rel[2],       0, -rel[0],
            -rel[1],  rel[0],       0);

    c = cscale * c;
    for(int jj=0; jj<2; jj++){
      for(int kk = 0; kk<2; kk++){
        K(row+jj, col+kk) = c(jj,kk);
        K(col+kk, row+jj) = c(jj,kk);
      }
    }*/ 

    for(int jj=0; jj<2; jj++){
      K(row+2+jj, col+jj) = 1;
      K(col+jj, row+2+jj) = 1;
    }
  }
}

int StepperNewton2D::compute_dx_dense(ElementMesh2D * iMesh, const std::vector<Vector2S> &iForces, bool iRmRigid, std::vector<cfgScalar> &bb)
{
  int ndof = (int)bb.size();
  assert(iMesh && 2*iMesh->x.size()==ndof);

  MatrixXS K = iMesh->getStiffness();
   
  for(unsigned int ii = 0;ii<m->x.size(); ii++){
    for(int jj = 0;jj<2;jj++){
      int row = 2*ii + jj;
      //damping, better condition number
  //    K(row,row) += 100;
      if(m->fixed[ii]){
        bb[ row ] = 0;
        for(int kk = 0;kk<ndof;kk++){
          K(row,kk) = 0;
          K(row,row) = 100;
        }
      }else{
        bb[ row ] = iForces[ii][jj];
      }
    }
  }
  if(iRmRigid){
    K.resize(ndof + 3, ndof+3);
    bb.resize(ndof+3);
    fixRigid(K,&(bb[0]),m);
  }
  //linSolvef(K,&(bb[0]));

  return 0;
}

int StepperNewton2D::compute_dx_sparse(ElementMesh2D * iMesh, const std::vector<Vector2S> &iForces, bool iRemoveTranslation, bool iRemoveRotation, bool iPeriodic, std::vector<cfgScalar> &bb)
{
  bool triangular = true;

  if (!m_Init)
  {
    m_I.clear();
    m_J.clear();
    iMesh->stiffnessPattern(m_I, m_J, triangular, iRemoveTranslation, iRemoveRotation, iPeriodic);
    sparseInit();
    m_Init = true;
  }

  assert(iMesh);
  int ndof = (int)(2*iMesh->x.size());
 
  std::vector<cfgScalar> Kvalues;
  iMesh->getStiffnessSparse(Kvalues, triangular, true, iRemoveTranslation, iRemoveRotation, iPeriodic);
   
  for(unsigned int ii = 0;ii<m->x.size(); ii++){
    for(int jj = 0;jj<2;jj++){
      int row = 2*ii + jj;
      if(m->fixed[ii]){
        bb[ row ] = 0;
      }else{
        bb[ row ] = iForces[ii][jj];
      }
    }
  }
  int irow, nrow=(int)bb.size();
  for (irow=ndof; irow<nrow; irow++)
  {
    bb [irow] = 0;
  }

  int ncol = nrow;
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

    //for (Eigen::SparseMatrix<cfgScalar>::InnerIterator it(K, col); it; ++it){
    //    int row = it.row();
    for (int i=m_I[col]; i<m_I[col+1]; i++)
    {
      int row = m_J[indVal];
      if(triangular && col>row)
      {
        continue;
      }
      ja.push_back(row);
      val.push_back(Kvalues[indVal]);
      indVal++;
    }
  }
  ia[nrow] = indVal;

  int status = sparseSolve( &(ia[0]), &(ja[0]), &(val[0]), nrow, &(x[0]), &(rhs[0]));
   //   std::cout<< "sparse solve " << status << "\n";
  for(int ii = 0;ii<x.size();ii++){
    bb[ii] = x[ii];
  }
  if (status<0)
    return 1;

  return 0;
}
