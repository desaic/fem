#include "ElementMesh.hpp"
#include "Element.hpp"
#include "Eigen/Sparse"
#include "Timer.hpp"
#include <iostream>

using namespace Eigen;

typedef Eigen::Triplet<cfgScalar> TripletS;

///@brief add rows to K and b to constrain 6 degrees of freedom.
void fixRigid(Eigen::SparseMatrix<cfgScalar> & K,  ElementMesh * mesh)
{
  int row = K.rows();
  K.conservativeResize(row + 6, K.cols()+6);
  Eigen::SparseMatrix<cfgScalar> KConstraint(row + 6,row + 6);

  Vector3S center = Vector3S::Zero();
  for(int ii = 0;ii<mesh->x.size();ii++){
    center += (Vector3S) mesh->x[ii];
  }
  center /= (cfgScalar)mesh->x.size();
  cfgScalar cScale = 1000;
  std::vector<TripletS> triplets;
  for(int ii = 0;ii<mesh->x.size();ii++){
    int col = 3*ii;
    //relative position to center
    Vector3S rel = (Vector3S) mesh->x[ii] - center;
    //cross product matrix
    Matrix3S c;
    c <<       0, -rel[2],  rel[1],
          rel[2],       0, -rel[0],
         -rel[1],  rel[0],       0;
    c = -cScale * c;
    for(int kk = 0; kk<3; kk++){
      triplets.push_back( TripletS(row + kk+3,   col + kk, -cScale ));
//      triplets.push_back( Tripf(col + kk  , row + kk+3, cScale ));
      for(int ll = 0; ll<3;ll++){
        triplets.push_back( TripletS(row + kk, col + ll, c(kk,ll)) );
//        triplets.push_back( Tripf(col + ll, row + kk, c(kk,ll)) );
      }
    }
  }
  for(int ii = 0;ii<6;ii++){
    triplets.push_back(TripletS(row+ii, row + ii,0.f));
  }
  KConstraint.setFromTriplets(triplets.begin(), triplets.end());
  K += KConstraint;
}

void ElementMesh::getStiffnessSparse(std::vector<cfgScalar> &val, bool trig, bool constrained, bool iFixedRigid)
{
  int N = 3* (int)x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N,N);
  for(unsigned int ii = 0;ii<e.size();ii++){
    Element * ele = e[ii];
    int nV = ele->nV();
    MatrixXS K  = getStiffness(ii);
    for(int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<3;dim1++){
          for(int dim2= 0 ;dim2<3;dim2++){
            if(trig && (3*vk+dim2 > 3*vj+dim1)) {
              continue;
            }
            cfgScalar val = K(3 * jj + dim1, 3 * kk + dim2);
            if (constrained){
              if ((fixed[vk] || fixed[vj]) && (3 * jj + dim1) != (3 * kk + dim2)){
                val = 0;
              }
            }
            TripletS triple(3 * vj + dim1, 3 * vk + dim2, val);
            coef.push_back(triple);
          }
        }
      }
    }
  }
  Ksparse.setFromTriplets(coef.begin(), coef.end());
  if (iFixedRigid)
  {
    fixRigid(Ksparse, this);
  }
  for(int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<cfgScalar>::InnerIterator it(Ksparse, ii); it; ++it){
     val.push_back(it.value());
   }
  }
 
}

Eigen::SparseMatrix<cfgScalar>
ElementMesh::getStiffnessSparse(bool trig, bool constrained, bool iFixedRigid)
{
  int N = 3* (int)x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N,N);
  for(unsigned int ii = 0;ii<e.size();ii++){
    Element * ele = e[ii];
    int nV = ele->nV();
    MatrixXS K  = getStiffness(ii);
    for(int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<3;dim1++){
          for(int dim2= 0 ;dim2<3;dim2++){
            if(trig && (3*vk+dim2 > 3*vj+dim1)) {
              continue;
            }
            cfgScalar val = K(3 * jj + dim1, 3 * kk + dim2);
            if (constrained){
              if ( (fixed[vk] || fixed[vj]) && (3 * jj + dim1 ) != (3 * kk + dim2) ){
                val = 0;
              }
            }
            TripletS triple(3 * vj + dim1, 3 * vk + dim2, val);
            coef.push_back(triple);
          }
        }
      }
    }
  }
  Ksparse.setFromTriplets(coef.begin(), coef.end());
  if (iFixedRigid)
  {
    fixRigid(Ksparse, this);
  }
  return Ksparse;
}

void ElementMesh::getStiffnessSparse(std::vector<cfgScalar> &val, bool trig, bool constrained, bool iFixedTranslation, bool iFixedRotation, bool iPeriodic)
{
  //Timer t;
  int N = 3* (int)x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N,N);
  //t.start();
  for(unsigned int ii = 0;ii<e.size();ii++){
    Element * ele = e[ii];
    int nV = ele->nV();
    MatrixXS K  = getStiffness(ii);
    for(int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<3;dim1++){
          for(int dim2= 0 ;dim2<3;dim2++){
            if(trig && (3*vk+dim2 > 3*vj+dim1)) {
              continue;
            }
            cfgScalar val = K(3 * jj + dim1, 3 * kk + dim2);
            if (constrained){
              if (fixed[vk] || fixed[vj]){
                val = 0;
                if (vj == vk && dim1 == dim2){
                  val = 100;
                }
              }
            }
            //if (vk == vj && dim1 == dim2){
            //  val += 1;
            //}
            TripletS triple(3 * vj + dim1, 3 * vk + dim2, val);
            coef.push_back(triple);
          }
        }
      }
    }
  }
  //t.end();
  //std::cout << " computeStiffness: " <<  t.getSeconds() << std::endl;

  //t.start();
  Ksparse.setFromTriplets(coef.begin(), coef.end());
  //t.end();
  //std::cout << " setFromTriplets: " <<  t.getSeconds() << std::endl;

  //t.start();
  if (iFixedTranslation && iFixedRotation)
  {
    fixTranslation(Ksparse, trig, this);
    fixRotation(Ksparse, trig, this);
  }
  else if (iFixedTranslation)
  {
    fixTranslation(Ksparse, trig, this);
  }
  else if (iFixedRotation)
  {
    fixRotation(Ksparse, trig, this);
  }
  if (iPeriodic)
  {
    enforcePeriodicity(Ksparse, trig, this);
  }

  for(int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<cfgScalar>::InnerIterator it(Ksparse, ii); it; ++it){
     val.push_back(it.value());
   }
  }
  //t.end();
  //std::cout << " fixConstraints + setvalues: " <<  t.getSeconds() << std::endl;
}

void ElementMesh::stiffnessPattern(std::vector<int> & I, std::vector<int> & J,
   bool trig, bool iFixedRigid)
{
  int N = 3* (int)x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N,N);

  for(unsigned int ii = 0;ii<e.size();ii++){
    Element * ele = e[ii];
    for(int jj = 0; jj<ele->nV(); jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<ele->nV(); kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<3;dim1++){
          for(int dim2= 0 ;dim2<3;dim2++){
            if(trig&&(3*vk+dim2 > 3*vj+dim1)){
              continue;
            }
            TripletS triple(3*vj+dim1,3*vk+dim2,1);
            coef.push_back(triple);
          }
        }
      }
    }
  }

  Ksparse.setFromTriplets(coef.begin(), coef.end());
  if (iFixedRigid)
  {
    fixRigid(Ksparse, this);
  }

  I.push_back(0);
  for(int ii = 0; ii<Ksparse.cols(); ii++){
    for (Eigen::SparseMatrix<cfgScalar>::InnerIterator it(Ksparse, ii); it; ++it){
     J.push_back(it.row());
   }
    I.push_back((int)J.size());
  }
}

void ElementMesh::stiffnessPattern(std::vector<int> & I, std::vector<int> & J,  bool trig, bool iFixedTranslation, bool iFixedRotation, bool iPeriodic)
{
  int N = 3 * (int)x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N,N);

  for(unsigned int ii = 0;ii<e.size();ii++){
    Element * ele = e[ii];
    for(int jj = 0; jj<ele->nV(); jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<ele->nV(); kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<3;dim1++){
          for(int dim2= 0 ;dim2<3;dim2++){
            if(trig&&(3*vk+dim2 > 3*vj+dim1)){
              continue;
            }
            TripletS triple(3*vj+dim1,3*vk+dim2,1);
            coef.push_back(triple);
          }
        }
      }
    }
  }

  Ksparse.setFromTriplets(coef.begin(), coef.end());
  if (iFixedTranslation && iFixedRotation)
  {
    fixTranslation(Ksparse, trig, this);
    fixRotation(Ksparse, trig, this);
  }
  else if (iFixedTranslation)
  {
    fixTranslation(Ksparse, trig, this);
  }
  else if (iFixedRotation)
  {
    fixRotation(Ksparse, trig, this);
  }
  if (iPeriodic)
  {
    enforcePeriodicity(Ksparse, trig, this);
  }

  I.push_back(0);
  for(int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<cfgScalar>::InnerIterator it(Ksparse, ii); it; ++it){
     J.push_back(it.row());
   }
    I.push_back((int)J.size());
  }
}

