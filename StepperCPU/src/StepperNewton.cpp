#include "StepperNewton.hpp"
#include "ElementMesh.hpp"
#include "MatrixX.hpp"
#include "ArrayUtil.hpp"
#include "femError.hpp"
#include "ElementRegGrid.hpp"

#include "LinSolve.hpp"
#include "SparseLin.hpp"
#include "Eigen/Sparse"

#include <fstream>
#include <assert.h>

StepperNewton::StepperNewton():dense(true),dx_tol(1e-5f),h(1.0f)
{
  m_Init = false;
  dense = false;

  m_periodic = false;
  m_noTranslation = false;
  m_noRotation = false;
}

void StepperNewton::enforcePeriodicity(bool iOn)
{
  m_periodic = iOn;
}

void StepperNewton::removeTranslation(bool iOn)
{
  m_noTranslation = iOn;
}

void StepperNewton::removeRotation(bool iOn)
{
  m_noRotation = iOn;
}

int StepperNewton::oneStep()
{
  int status = 0;

  std::vector<Vector3f> force = m->getForce();
  float E = m->getEnergy();

  int ndof = 3*(int)m->x.size();
  int nconstraint = 0;
  if (m_noTranslation)
  {
    nconstraint +=3;
  }
  if (m_noRotation)
  {
    nconstraint +=3;
  }
  if (m_periodic)
  {
    int nx = ((ElementRegGrid*)m)->nx+1;
    int ny = ((ElementRegGrid*)m)->ny+1;
    int nz = ((ElementRegGrid*)m)->nz+1;
    nconstraint += 3*(ny*nz-1 + (nx-1)*nz-1 + (nx-1)*(ny-1)-1);
  }
  std::vector<float> bb(ndof+nconstraint);

  bool rmRigid = m_noTranslation && m_noRotation;

  if (dense)
  {
    status = compute_dx_dense(m, force, rmRigid, bb);
  }
  else
  {
    if (!m_periodic)
    {
      status = compute_dx_sparse(m, force, rmRigid, bb);
    }
    else
    {
      status = compute_dx_sparse(m, force, m_noTranslation, m_noRotation, m_periodic, bb);
    }
  }
  if (status>0)
    return status;

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
     // h=1.1f*h;
      break;
    }
  }
  if (std::abs(E - E1) < 1e-3f){
    return -1;
  }
  std::cout << E << "\n";
  return 0;
}

///@brief add rows to K and b to constrain 6 degrees of freedom.
///@param K size is #DOF + 6
void fixRigid(MatrixXf & K, float * b,
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

int StepperNewton::compute_dx_dense(ElementMesh * iMesh, const std::vector<Vector3f> &iForces, bool iRmRigid, std::vector<float> &bb)
{
  int ndof = bb.size();
  assert(iMesh && 3*iMesh->x.size()==ndof);

  MatrixXf K = iMesh->getStiffness();
   
  for(unsigned int ii = 0;ii<m->x.size(); ii++){
    for(int jj = 0;jj<3;jj++){
      int row = 3*ii + jj;
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
    K.resize(ndof + 6, ndof+6);
    bb.resize(ndof+6);
    fixRigid(K,&(bb[0]),m);
  }
  linSolvef(K,&(bb[0]));

  return 0;
}

int StepperNewton::compute_dx_sparse(ElementMesh * iMesh, const std::vector<Vector3f> &iForces, bool iRmRigid, std::vector<float> &bb)
{
  bool triangular = true;

  if (!m_Init)
  {
    m_I.clear();
    m_J.clear();
    iMesh->stiffnessPattern(m_I, m_J, triangular, iRmRigid);
    sparseInit();
    m_Init = true;
  }

  int ndof = bb.size();
  assert(iMesh && 3*iMesh->x.size()==bb.size());

  std::vector<float> Kvalues;
  iMesh->getStiffnessSparse(Kvalues, triangular, true, iRmRigid);
   
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
  if (status<0)
    return 1;

  return 0;
}

int StepperNewton::compute_dx_sparse(ElementMesh * iMesh, const std::vector<Vector3f> &iForces, bool iRemoveTranslation, bool iRemoveRotation, bool iPeriodic, std::vector<float> &bb)
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
  int ndof = (int)(3*iMesh->x.size());
 
  std::vector<float> Kvalues;
  iMesh->getStiffnessSparse(Kvalues, triangular, true, iRemoveTranslation, iRemoveRotation, iPeriodic);
   
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
