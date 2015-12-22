#include "ElementMesh.hpp"
#include "Element.hpp"
#include "Eigen/Sparse"

using namespace Eigen;

typedef Eigen::Triplet<float> Tripletf;

///@brief add rows to K and b to constrain 6 degrees of freedom.
void fixRigid(Eigen::SparseMatrix<float> & K,  ElementMesh * mesh)
{
  int row = K.rows();
  K.conservativeResize(row + 6, K.cols()+6);
  Eigen::SparseMatrix<float> KConstraint(row + 6,row + 6);

  Eigen::Vector3f center = Eigen::Vector3f::Zero();
  for(int ii = 0;ii<mesh->x.size();ii++){
    center += (Eigen::Vector3f) mesh->x[ii];
  }
  center /= (float)mesh->x.size();
  float cScale = 1000;
  std::vector<Tripletf> triplets;
  for(int ii = 0;ii<mesh->x.size();ii++){
    int col = 3*ii;
    //relative position to center
    Eigen::Vector3f rel = (Eigen::Vector3f) mesh->x[ii] - center;
    //cross product matrix
    Eigen::Matrix3f c;
    c <<       0, -rel[2],  rel[1],
          rel[2],       0, -rel[0],
         -rel[1],  rel[0],       0;
    c = -cScale * c;
    for(int kk = 0; kk<3; kk++){
      triplets.push_back( Tripletf(row + kk+3,   col + kk, -cScale ));
//      triplets.push_back( Tripf(col + kk  , row + kk+3, cScale ));
      for(int ll = 0; ll<3;ll++){
        triplets.push_back( Tripletf(row + kk, col + ll, c(kk,ll)) );
//        triplets.push_back( Tripf(col + ll, row + kk, c(kk,ll)) );
      }
    }
  }
  for(int ii = 0;ii<6;ii++){
    triplets.push_back(Tripletf(row+ii, row + ii,0.f));
  }
  KConstraint.setFromTriplets(triplets.begin(), triplets.end());
  K += KConstraint;
}

void ElementMesh::getStiffnessSparse(std::vector<float> &val, bool trig, bool constrained, bool iFixedRigid)
{
  int N = 3* (int)x.size();
  std::vector<Tripletf> coef;
  Eigen::SparseMatrix<float> Ksparse(N,N);
  for(unsigned int ii = 0;ii<e.size();ii++){
    Element * ele = e[ii];
    int nV = ele->nV();
    MatrixXf K  = getStiffness(ii);
    for(int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<3;dim1++){
          for(int dim2= 0 ;dim2<3;dim2++){
            if(trig && (3*vk+dim2 > 3*vj+dim1)) {
              continue;
            }
            float val = K(3 * jj + dim1, 3 * kk + dim2);
            if (constrained){
              if ((fixed[vk] || fixed[vj]) && (3 * jj + dim1) != (3 * kk + dim2)){
                val = 0;
              }
            }
            Tripletf triple(3 * vj + dim1, 3 * vk + dim2, val);
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
    for (Eigen::SparseMatrix<float>::InnerIterator it(Ksparse, ii); it; ++it){
     val.push_back(it.value());
   }
  }
 
}

Eigen::SparseMatrix<float>
ElementMesh::getStiffnessSparse(bool trig, bool constrained, bool iFixedRigid)
{
  int N = 3* (int)x.size();
  std::vector<Tripletf> coef;
  Eigen::SparseMatrix<float> Ksparse(N,N);
  for(unsigned int ii = 0;ii<e.size();ii++){
    Element * ele = e[ii];
    int nV = ele->nV();
    MatrixXf K  = getStiffness(ii);
    for(int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<3;dim1++){
          for(int dim2= 0 ;dim2<3;dim2++){
            if(trig && (3*vk+dim2 > 3*vj+dim1)) {
              continue;
            }
            float val = K(3 * jj + dim1, 3 * kk + dim2);
            if (constrained){
              if ( (fixed[vk] || fixed[vj]) && (3 * jj + dim1 ) != (3 * kk + dim2) ){
                val = 0;
              }
            }
            Tripletf triple(3 * vj + dim1, 3 * vk + dim2, val);
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

void ElementMesh::getStiffnessSparse(std::vector<float> &val, bool trig, bool constrained, bool iFixedTranslation, bool iFixedRotation, bool iPeriodic)
{
  int N = 3* (int)x.size();
  std::vector<Tripletf> coef;
  Eigen::SparseMatrix<float> Ksparse(N,N);
  for(unsigned int ii = 0;ii<e.size();ii++){
    Element * ele = e[ii];
    int nV = ele->nV();
    MatrixXf K  = getStiffness(ii);
    for(int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<3;dim1++){
          for(int dim2= 0 ;dim2<3;dim2++){
            if(trig && (3*vk+dim2 > 3*vj+dim1)) {
              continue;
            }
            float val = K(3 * jj + dim1, 3 * kk + dim2);
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
            Tripletf triple(3 * vj + dim1, 3 * vk + dim2, val);
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

  for(int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<float>::InnerIterator it(Ksparse, ii); it; ++it){
     val.push_back(it.value());
   }
  }
}

void ElementMesh::stiffnessPattern(std::vector<int> & I, std::vector<int> & J,
   bool trig, bool iFixedRigid)
{
  int N = 3* (int)x.size();
  std::vector<Tripletf> coef;
  Eigen::SparseMatrix<float> Ksparse(N,N);

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
            Tripletf triple(3*vj+dim1,3*vk+dim2,1);
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
    for (Eigen::SparseMatrix<float>::InnerIterator it(Ksparse, ii); it; ++it){
     J.push_back(it.row());
   }
    I.push_back((int)J.size());
  }
}

void ElementMesh::stiffnessPattern(std::vector<int> & I, std::vector<int> & J,  bool trig, bool iFixedTranslation, bool iFixedRotation, bool iPeriodic)
{
  int N = 3 * (int)x.size();
  std::vector<Tripletf> coef;
  Eigen::SparseMatrix<float> Ksparse(N,N);

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
            Tripletf triple(3*vj+dim1,3*vk+dim2,1);
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
    for (Eigen::SparseMatrix<float>::InnerIterator it(Ksparse, ii); it; ++it){
     J.push_back(it.row());
   }
    I.push_back((int)J.size());
  }
}

