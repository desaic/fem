#include "ElementMesh2D.h"
#include "Material2D.h"
#include "Element2D.h"
#include "femError.hpp"
#include <iostream>
#include "MeshUtilities.h"
#include "ElementRegGrid2D.h"

typedef Eigen::Triplet<cfgScalar> TripletS;

void ElementMesh2D::initArrays()
{
  x=X;
  fe.resize(X.size(), Vector2S(0,0));
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

cfgScalar ElementMesh2D::getEnergy()
{
  cfgScalar ene = 0;
  for(unsigned int ii = 0;ii<e.size();ii++){
    ene += getEnergy(ii);
    if(fem_error){
      return -1;
    }
  }
  //energy from external forces
  for(unsigned int ii = 0;ii<fe.size();ii++){
    ene -= fe[ii].dot(x[ii]);
  }
  return ene;
}

cfgScalar ElementMesh2D::getEnergy(int eIdx)
{
  return m[me[eIdx]]->getEnergy(e[eIdx], this);
}

std::vector<Vector2S> ElementMesh2D::getForce(int eIdx)
{
  return m[me[eIdx]]->getForce(e[eIdx], this);
}

std::vector<Vector2S> ElementMesh2D::getForce()
{
  std::vector<Vector2S> force(x.size(), Vector2S(0,0));
  for(unsigned int ii = 0;ii<e.size();ii++){
    std::vector<Vector2S> fele = getForce(ii);
    for(int jj = 0; jj<e[ii]->nV(); jj++){
      force[ e[ii]->at(jj) ] += fele[jj];
    }
  }
  for(unsigned int ii= 0;ii<fe.size();ii++){
    force[ii] += fe[ii];
  }
  return force;
}

cfgScalar ElementMesh2D::eleSize()
{
  return X[e[0]->at(7)][0] - X[e[0]->at(0)][0];
}

ElementMesh2D::ElementMesh2D() :u(0), forceDrawingScale(10)
{}

ElementMesh2D::~ElementMesh2D()
{
  for(unsigned int ii = 0; ii<e.size();ii++){
    delete e[ii];
  }
}


///@brief add rows to K and b to constrain 6 degrees of freedom.
void ElementMesh2D::fixRigid(Eigen::SparseMatrix<cfgScalar> & K,  bool iTriangular, ElementMesh2D * mesh)
{
  int nrow = K.rows();
  int ncol = K.cols();
  int nconstraints = 3;

  K.conservativeResize(nrow + nconstraints, ncol + nconstraints);
  Eigen::SparseMatrix<cfgScalar> KConstraint(nrow + nconstraints, nrow + nconstraints);

  Vector2S center = Vector2S::Zero();
  for(int ii = 0;ii<mesh->x.size();ii++){
    center += (Vector2S) mesh->x[ii];
  }
  center /= (cfgScalar)mesh->x.size();
  //cfgScalar cScale = 1000;
  cfgScalar cScale = 1;
  std::vector<TripletS> triplets;
  for(int ii = 0;ii<mesh->x.size();ii++)
  {
    int col = 2*ii;
    //relative position to center
    Vector2S rel = (Vector2S) mesh->x[ii] - center;
    Vector2S c(-rel[1], rel[0]);  // no rotation
    c = -cScale * c;

     //fix translation
    for(int kk=0; kk<2; kk++)
    {
      assert(!iTriangular || nrow+kk >= col+kk);
      triplets.push_back( TripletS(nrow+kk,   col+kk, -cScale));    
      if (!iTriangular)
      {
        triplets.push_back( TripletS(col+kk, nrow+kk, cScale));    
      }
    }
    int shift = 2;

    // fix rotation
    for(int ll=0; ll<2;ll++)
    {
      assert(!iTriangular || nrow+shift >= col + ll);
      triplets.push_back( TripletS(nrow+shift, col + ll, c(ll)) );
      if (!iTriangular)
      {
        triplets.push_back( TripletS(col + ll, nrow+shift, -c(ll)));
      }
    }
  }
  for(int ii = 0;ii<nconstraints;ii++){
    triplets.push_back(TripletS(nrow+ii, nrow+ii,0.f));
  }
  KConstraint.setFromTriplets(triplets.begin(), triplets.end());
  K += KConstraint;
}

void ElementMesh2D::fixTranslation(Eigen::SparseMatrix<cfgScalar> & K,  bool iTriangular, ElementMesh2D * mesh)
{
  int nrow = K.rows();
  int ncol = K.cols();
  int nconstraints = 2;

  K.conservativeResize(nrow + nconstraints, ncol + nconstraints);
  Eigen::SparseMatrix<cfgScalar> KConstraint(nrow + nconstraints, nrow + nconstraints);

  Vector2S center = Vector2S::Zero();
  for(int ii = 0;ii<mesh->x.size();ii++){
    center += (Vector2S) mesh->x[ii];
  }
  center /= (cfgScalar)mesh->x.size();
  cfgScalar cScale = 1;
  std::vector<TripletS> triplets;
  for(int ii = 0;ii<mesh->x.size();ii++)
  {
    int col = 2*ii;
    //relative position to center
    Vector2S rel = (Vector2S) mesh->x[ii] - center;
   
     //fix translation
    for(int kk=0; kk<2; kk++)
    {
      assert(!iTriangular || nrow+kk >= col+kk);
      triplets.push_back( TripletS(nrow+kk,   col+kk, -cScale));    
      if (!iTriangular)
      {
        triplets.push_back( TripletS(col+kk, nrow+kk, cScale));    
      }
    }
  }
  for(int ii = 0;ii<nconstraints;ii++){
    triplets.push_back(TripletS(nrow+ii, nrow+ii,0.f));
  }
  KConstraint.setFromTriplets(triplets.begin(), triplets.end());
  K += KConstraint;
}

void ElementMesh2D::fixRotation(Eigen::SparseMatrix<cfgScalar> & K,  bool iTriangular, ElementMesh2D * mesh)
{
  int nrow = K.rows();
  int ncol = K.cols();
  int nconstraints = 1;

  K.conservativeResize(nrow + nconstraints, ncol + nconstraints);
  Eigen::SparseMatrix<cfgScalar> KConstraint(nrow + nconstraints, nrow + nconstraints);

  Vector2S center = Vector2S::Zero();
  for(int ii = 0;ii<mesh->x.size();ii++){
    center += (Vector2S) mesh->x[ii];
  }
  center /= (cfgScalar)mesh->x.size();
  cfgScalar cScale = 1;
  std::vector<TripletS> triplets;
  for(int ii = 0;ii<mesh->x.size();ii++)
  {
    int col = 2*ii;
    //relative position to center
    Vector2S rel = (Vector2S) mesh->x[ii] - center;
    Vector2S c(-rel[1], rel[0]);  // no rotation
    c = -cScale * c;

    // fix rotation
    for(int ll=0; ll<2;ll++)
    {
      assert(!iTriangular || nrow >= col + ll);
      triplets.push_back( TripletS(nrow, col + ll, c(ll)) );
      if (!iTriangular)
      {
        triplets.push_back( TripletS(col + ll, nrow, -c(ll)));
      }
    }
  }
  for(int ii = 0;ii<nconstraints;ii++){
    triplets.push_back(TripletS(nrow+ii, nrow+ii,0.f));
  }
  KConstraint.setFromTriplets(triplets.begin(), triplets.end());
  K += KConstraint;
}

void ElementMesh2D::enforcePeriodicity(Eigen::SparseMatrix<cfgScalar> & K, bool iTriangular, ElementMesh2D * mesh)
{
  std::vector<int> sideVertexIndices[4];
  int iside;
  for (iside=0; iside<4; iside++)
  {
    meshUtil::getSideVertices(iside, (const ElementRegGrid2D*)this, sideVertexIndices[iside]);
  }
  int nconstraints = 2*((int)sideVertexIndices[1].size()-1 + (int)sideVertexIndices[3].size()-2);

  int nrow = K.rows();
  int ncol = K.cols();

  int nrowInit = nrow;
  int ncolInit = ncol;

  K.conservativeResize(nrow + nconstraints, ncol + nconstraints);
  Eigen::SparseMatrix<cfgScalar> KConstraint(nrow + nconstraints, nrow + nconstraints);

  std::vector<TripletS> triplets;
  cfgScalar cScale = 1;

  int indCorner0 = sideVertexIndices[0][0];
  int indCorner1 = sideVertexIndices[1][0]; // right
  int indCorner3 = sideVertexIndices[3][0]; // top

  int ivertex, nvertex=(int)sideVertexIndices[1].size();
  for (ivertex=1; ivertex<nvertex; ivertex++)
  {
    int indVertex0 =  sideVertexIndices[0][ivertex];
    int indVertex1 =  sideVertexIndices[1][ivertex];

    for (int icoord=0; icoord<2; icoord++)
    {
      int row = nrow+icoord;

      int col0 = 2*indCorner0 + icoord;
      int col1 = 2*indCorner1 + icoord;
      int col2 = 2*indVertex0 + icoord;
      int col3 = 2*indVertex1 + icoord;

      assert(!iTriangular || (row >= col0 && row >= col1 && row>= col2 && row>= col3));
      triplets.push_back( TripletS(row,   col0, -cScale)); 
      triplets.push_back( TripletS(row,   col1, cScale)); 
      triplets.push_back( TripletS(row,   col2, cScale)); 
      triplets.push_back( TripletS(row,   col3, -cScale)); 
      if (!iTriangular)
      {
        triplets.push_back( TripletS(col0, row, cScale));
        triplets.push_back( TripletS(col1, row, -cScale));    
        triplets.push_back( TripletS(col2, row, -cScale));    
        triplets.push_back( TripletS(col3, row, cScale));    
      }
    }
    nrow += 2;
  }
  nvertex = (int)sideVertexIndices[3].size()-1;
  for (ivertex=1; ivertex<nvertex; ivertex++)
  {
    int indVertex0 =  sideVertexIndices[2][ivertex];
    int indVertex1 =  sideVertexIndices[3][ivertex];

    for (int icoord=0; icoord<2; icoord++)
    {
      int row = nrow+icoord;

      int col0 = 2*indCorner0 + icoord;
      int col1 = 2*indCorner3 + icoord;
      int col2 = 2*indVertex0 + icoord;
      int col3 = 2*indVertex1 + icoord;

      assert(!iTriangular || (row >= col0 && row >= col1 && row>= col2 && row>= col3));
      triplets.push_back( TripletS(row,   col0, -cScale)); 
      triplets.push_back( TripletS(row,   col1, cScale)); 
      triplets.push_back( TripletS(row,   col2, cScale)); 
      triplets.push_back( TripletS(row,   col3, -cScale)); 
      if (!iTriangular)
      {
        triplets.push_back( TripletS(col0, row, cScale));
        triplets.push_back( TripletS(col1, row, -cScale));    
        triplets.push_back( TripletS(col2, row, -cScale));    
        triplets.push_back( TripletS(col3, row, cScale));    
      }
    }
    nrow += 2;
  }

  for(int ii = 0;ii<nconstraints;ii++){
    triplets.push_back(TripletS(nrowInit+ii, nrowInit+ii,0.f));
  }
  KConstraint.setFromTriplets(triplets.begin(), triplets.end());
  K += KConstraint;
}

void ElementMesh2D::getStiffnessSparse(std::vector<cfgScalar> &val, bool trig, bool constrained, bool iFixedTranslation, bool iFixedRotation, bool iPeriodic, Eigen::SparseMatrix<cfgScalar> *oMatrix)
{
  int N = 2* (int)x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N,N);
  for(unsigned int ii = 0;ii<e.size();ii++){
    Element2D * ele = e[ii];
    int nV = ele->nV();
    MatrixXS K  = getStiffness(ii);
    for(int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for(int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for(int dim1= 0 ;dim1<2;dim1++){
          for(int dim2= 0 ;dim2<2;dim2++){
            if(trig && (2*vk+dim2 > 2*vj+dim1)) {
              continue;
            }
            cfgScalar val = K(2 * jj + dim1, 2 * kk + dim2);
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
            TripletS triple(2 * vj + dim1, 2 * vk + dim2, val);
            coef.push_back(triple);
          }
        }
      }
    }
  }
  Ksparse.setFromTriplets(coef.begin(), coef.end());
  //if (iFixedRigid)
  if (iFixedTranslation && iFixedRotation)
  {
    //fixRigid(Ksparse, trig, this);
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
  if (oMatrix)
  {
    *oMatrix = Ksparse;
  }

  for(int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<cfgScalar>::InnerIterator it(Ksparse, ii); it; ++it){
     val.push_back(it.value());
   }
  }
 
}

void ElementMesh2D::stiffnessPattern(std::vector<int> & I, std::vector<int> & J,  bool trig, bool iFixedTranslation, bool iFixedRotation, bool iPeriodic)
{
  int N = 2* (int)x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N,N);

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
            TripletS triple(2*vj+dim1,2*vk+dim2,1);
            coef.push_back(triple);
          }
        }
      }
    }
  }

  Ksparse.setFromTriplets(coef.begin(), coef.end());
   //if (iFixedRigid)
  if (iFixedTranslation && iFixedRotation)
  {
    //fixRigid(Ksparse, trig, this);
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


MatrixXS ElementMesh2D::getStiffness(int eIdx)
{
  MatrixXS K = m[me[eIdx]]->getStiffness(e[eIdx],this);
  return K;
}

MatrixXS ElementMesh2D::getStiffness()
{
  int matSize = 2 * (int)x.size();
  MatrixXS Kglobal=MatrixXS::Zero(matSize,matSize);
  for(unsigned int ii = 0;ii<e.size();ii++){
    MatrixXS K = getStiffness(ii);
    for(int jj = 0; jj<e[ii]->nV(); jj++){
      int vj = e[ii]->at(jj);
      for(int kk = 0; kk<e[ii]->nV(); kk++){
        int vk = e[ii]->at(kk);
        Kglobal.block(2 * vj, 2 * vk, 2, 2) = K.block(2 * jj, 2 * kk, 2, 2);
      }
    }
  }
  return Kglobal;
}

