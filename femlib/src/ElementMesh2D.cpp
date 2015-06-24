#include "ElementMesh2D.h"
#include "Material2D.h"
#include "Element2D.h"
#include "femError.hpp"
#include <iostream>

typedef Eigen::Triplet<float> Tripletf;

void ElementMesh2D::initArrays()
{
  x=X;
  fe.resize(X.size());
  fixed.resize(X.size());
  me.resize(e.size());
}

void ElementMesh2D::addMaterial(Material2D*_m)
{
  m.push_back(_m);
  _m->init(this);
}

int ElementMesh2D::check()
{
  if(fe.size()!=x.size()){
    std::cout<<"external force and dof size differ\n";
    return -1;
  }
  if(e.size()!=me.size()){
    std::cout<<"material assignment and num elements differ\n";
    return -1;
  }
  for(unsigned int ii = 0;ii<me.size();ii++){
    if(me[ii]<0 || me[ii]>=m.size()){
      std::cout<<"Material index out of bounds\n";
      std::cout<<"ele: "<<ii<<", mat: "<<me[ii]<<"\n";
      return -1;
    }
  }
    
  if(fixed.size() != x.size()){
    std::cout<<"fixed array differ in size to vertex array \n";
    return -1;
  }
  return 0;
}

float ElementMesh2D::getEnergy()
{
  float ene = 0;
  for(unsigned int ii = 0;ii<e.size();ii++){
    ene += getEnergy(ii);
    if(fem_error){
      return -1;
    }
  }
  //energy from external forces
  for(unsigned int ii = 0;ii<fe.size();ii++){
    ene -= Vector2f::dot(fe[ii], x[ii]);
  }
  return ene;
}

float ElementMesh2D::getEnergy(int eIdx)
{
  return m[me[eIdx]]->getEnergy(e[eIdx], this);
}

std::vector<Vector2f> ElementMesh2D::getForce(int eIdx)
{
  return m[me[eIdx]]->getForce(e[eIdx], this);
}

std::vector<Vector2f> ElementMesh2D::getForce()
{
  std::vector<Vector2f> force(x.size());
  for(unsigned int ii = 0;ii<e.size();ii++){
    std::vector<Vector2f> fele = getForce(ii);
    for(int jj = 0; jj<e[ii]->nV(); jj++){
      force[ e[ii]->at(jj) ] += fele[jj];
    }
  }
  for(unsigned int ii= 0;ii<fe.size();ii++){
    force[ii] += fe[ii];
  }
  return force;
}

float ElementMesh2D::eleSize()
{
  return X[e[0]->at(7)][0] - X[e[0]->at(0)][0];
}

ElementMesh2D::ElementMesh2D():u(0)
{}

ElementMesh2D::~ElementMesh2D()
{
  for(unsigned int ii = 0; ii<e.size();ii++){
    delete e[ii];
  }
}


///@brief add rows to K and b to constrain 6 degrees of freedom.
void ElementMesh2D::fixRigid(Eigen::SparseMatrix<float> & K,  ElementMesh2D * mesh)
{
  int row = K.rows();
  K.conservativeResize(row + 3, K.cols()+3);
  Eigen::SparseMatrix<float> KConstraint(row + 3,row + 3);

  Eigen::Vector2f center = Eigen::Vector2f::Zero();
  for(int ii = 0;ii<mesh->x.size();ii++){
    center += (Eigen::Vector2f) mesh->x[ii];
  }
  center /= mesh->x.size();
  float cScale = 1000;
  std::vector<Tripletf> triplets;
  for(int ii = 0;ii<mesh->x.size();ii++){
    int col = 2*ii;
    //relative position to center
    Eigen::Vector2f rel = (Eigen::Vector2f) mesh->x[ii] - center;
    //cross product matrix
   /* Eigen::Matrix3f c;
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
    }*/ 
  }
  for(int ii = 0;ii<3;ii++){
    triplets.push_back(Tripletf(row+ii, row + ii,0.f));
  }
  KConstraint.setFromTriplets(triplets.begin(), triplets.end());
  K += KConstraint;
}

void ElementMesh2D::getStiffnessSparse(std::vector<float> &val, bool trig, bool constrained, bool iFixedRigid)
{
  int N = 2* (int)x.size();
  std::vector<Tripletf> coef;
  Eigen::SparseMatrix<float> Ksparse(N,N);
  for(unsigned int ii = 0;ii<e.size();ii++){
    Element2D * ele = e[ii];
    int nV = ele->nV();
    MatrixXf K  = getStiffness(ii);
    for(int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<2;dim1++){
          for(int dim2= 0 ;dim2<2;dim2++){
            if(trig && (2*vk+dim2 > 2*vj+dim1)) {
              continue;
            }
            float val = K(2 * jj + dim1, 2 * kk + dim2);
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
            Tripletf triple(2 * vj + dim1, 2 * vk + dim2, val);
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

void ElementMesh2D::stiffnessPattern(std::vector<int> & I, std::vector<int> & J,  bool trig, bool iFixedRigid)
{
  int N = 2* (int)x.size();
  std::vector<Tripletf> coef;
  Eigen::SparseMatrix<float> Ksparse(N,N);

  for(unsigned int ii = 0;ii<e.size();ii++){
    Element2D * ele = e[ii];
    for(int jj = 0; jj<ele->nV(); jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<ele->nV(); kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<2;dim1++){
          for(int dim2= 0 ;dim2<2;dim2++){
            if(trig&&(2*vk+dim2 > 2*vj+dim1)){
              continue;
            }
            Tripletf triple(2*vj+dim1,2*vk+dim2,1);
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
  for(int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<float>::InnerIterator it(Ksparse, ii); it; ++it){
     J.push_back(it.row());
   }
    I.push_back((int)J.size());
  }
}


MatrixXf ElementMesh2D::getStiffness(int eIdx)
{
  MatrixXf K = m[me[eIdx]]->getStiffness(e[eIdx],this);
  return K;
}

MatrixXf ElementMesh2D::getStiffness()
{
  int matSize = 2 * (int)x.size();
  MatrixXf Kglobal(matSize,matSize);
  Kglobal.fill(0);
  for(unsigned int ii = 0;ii<e.size();ii++){
    MatrixXf K = getStiffness(ii);
    for(int jj = 0; jj<e[ii]->nV(); jj++){
      int vj = e[ii]->at(jj);
      for(int kk = 0; kk<e[ii]->nV(); kk++){
        int vk = e[ii]->at(kk);    
        addSubMat(K,Kglobal, 2*jj, 2*kk, 2*vj, 2*vk, 2, 2);
      }
    }
  }
  return Kglobal;
}

