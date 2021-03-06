#include "Element.hpp"
#include "ElementMesh.hpp"
#include "Material.hpp"
#include "ElementRegGrid.hpp"
#include "MeshUtilities.h"

using namespace Eigen;

typedef Eigen::Triplet<cfgScalar> TripletS;

MatrixXS ElementMesh::getStiffness(int eIdx)
{
  MatrixXS K = m[me[eIdx]]->getStiffness(e[eIdx],this);
  return K;
}

MatrixXS ElementMesh::getStiffness()
{
  int matSize = 3 * (int)x.size();
  MatrixXS Kglobal = MatrixXS::Zero(matSize,matSize);
  for(unsigned int ii = 0;ii<e.size();ii++){
    MatrixXS K = getStiffness(ii);
    for(int jj = 0; jj<e[ii]->nV(); jj++){
      int vj = e[ii]->at(jj);
      for(int kk = 0; kk<e[ii]->nV(); kk++){
        int vk = e[ii]->at(kk);    
        Kglobal.block(3 * vj, 3 * vk, 3, 3) = K.block(3 * jj, 3 * kk, 3, 3);
      }
    }
  }
  return Kglobal;
}


void getEleX(int ii, const ElementMesh * m, std::vector<Vector3S> &x)
{
  Element * ele = m->e[ii];
  x.resize(ele->nV());
  for (int jj = 0; jj<ele->nV(); jj++){
    x[jj] = m->x[ele->at(jj)];
  }
}

void setEleX(int ii, ElementMesh * m, const std::vector<Vector3S> &x)
{
  Element * ele = m->e[ii];
  for (int jj = 0; jj<ele->nV(); jj++){
    m->x[ele->at(jj)] = x[jj];
  }
}


void ElementMesh::fixTranslation(Eigen::SparseMatrix<cfgScalar> & K,  bool iTriangular, ElementMesh * mesh)
{
  int nrow = K.rows();
  int ncol = K.cols();
  int nconstraints = 3;

  K.conservativeResize(nrow + nconstraints, ncol + nconstraints);
  Eigen::SparseMatrix<cfgScalar> KConstraint(nrow + nconstraints, nrow + nconstraints);

  Vector3S center(0,0,0);
  for(int ii = 0;ii<mesh->x.size();ii++){
    center += mesh->x[ii];
  }
  center /= (cfgScalar)mesh->x.size();
  cfgScalar cScale = 1;
  std::vector<TripletS> triplets;
  for(int ii = 0;ii<mesh->x.size();ii++)
  {
    int col = 3*ii;
    //relative position to center
    Vector3S rel = mesh->x[ii] - center;
   
     //fix translation
    for(int kk=0; kk<3; kk++)
    {
      assert(!iTriangular || nrow+kk >= col+kk);
      triplets.push_back( TripletS(nrow+kk,   col+kk, -cScale));    
      if (!iTriangular)
      {
        triplets.push_back(TripletS(col+kk, nrow+kk, cScale));    
      }
    }
  }
  for(int ii = 0;ii<nconstraints;ii++){
    triplets.push_back(TripletS(nrow+ii, nrow+ii,0.f));
  }
  KConstraint.setFromTriplets(triplets.begin(), triplets.end());
  K += KConstraint;
}

void ElementMesh::fixRotation(Eigen::SparseMatrix<cfgScalar> & K,  bool iTriangular, ElementMesh * mesh)
{
  int nrow = K.rows();
  int ncol = K.cols();
  int nconstraints = 3;

  K.conservativeResize(nrow + nconstraints, ncol + nconstraints);
  Eigen::SparseMatrix<cfgScalar> KConstraint(nrow + nconstraints, nrow + nconstraints);

  Vector3S center(0,0,0);
  for(int ii = 0;ii<mesh->x.size();ii++){
    center += mesh->x[ii];
  }
  center /= mesh->x.size();
  cfgScalar cScale = 1;
  std::vector<TripletS> triplets;
  for(int ii = 0;ii<mesh->x.size();ii++)
  {
    int col = 3*ii;
    //relative position to center
    Vector3S rel = mesh->x[ii] - center;
    Matrix3S c;
    c <<       0, -rel[2],  rel[1],
          rel[2],       0, -rel[0],
         -rel[1],  rel[0],       0;
    c = -cScale * c;

    // fix rotation
    for(int kk = 0; kk<3; kk++)
    {
      for(int ll = 0; ll<3;ll++)
      {
        assert(!iTriangular || nrow >= col + ll);
        triplets.push_back( TripletS(nrow + kk, col + ll, c(kk,ll)) );
        if (!iTriangular)
        {
          triplets.push_back( TripletS(col + ll, nrow + kk, -c(kk,ll)));
        }
      }
    }

  }
  for(int ii = 0;ii<nconstraints;ii++){
    triplets.push_back(TripletS(nrow+ii, nrow+ii,0.f));
  }
  KConstraint.setFromTriplets(triplets.begin(), triplets.end());
  K += KConstraint;
}

void ElementMesh::enforcePeriodicity(Eigen::SparseMatrix<cfgScalar> & K, bool iTriangular, ElementMesh * mesh)
{
  std::vector<int> sideVertexIndices[6];
  int iside;
  for (iside=0; iside<6; iside++)
  {
    meshUtil::getSideVertices(iside, (const ElementRegGrid*)this, sideVertexIndices[iside]);
  }
  int nx = ((ElementRegGrid*)this)->nx+1;
  int ny = ((ElementRegGrid*)this)->ny+1;
  int nz = ((ElementRegGrid*)this)->nz+1;
  int nconstraints = 3*(ny*nz-1 + (nx-1)*nz-1 + (nx-1)*(ny-1)-1);

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
  int indCorner5 = sideVertexIndices[5][0]; // front

  int ivertex, nvertex=(int)sideVertexIndices[1].size();
  for (ivertex=1; ivertex<nvertex; ivertex++)
  {
    int indVertex0 =  sideVertexIndices[0][ivertex];
    int indVertex1 =  sideVertexIndices[1][ivertex];

    for (int icoord=0; icoord<3; icoord++)
    {
      int row = nrow+icoord;

      int col0 = 3*indCorner0 + icoord;
      int col1 = 3*indCorner1 + icoord;
      int col2 = 3*indVertex0 + icoord;
      int col3 = 3*indVertex1 + icoord;

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
    nrow += 3;
  }
  nvertex = (int)sideVertexIndices[3].size();
  for (ivertex=1; ivertex<nvertex; ivertex++)
  {
    if (ivertex%nx==nx-1) 
      continue;

    int indVertex0 =  sideVertexIndices[2][ivertex];
    int indVertex1 =  sideVertexIndices[3][ivertex];

    for (int icoord=0; icoord<3; icoord++)
    {
      int row = nrow+icoord;

      int col0 = 3*indCorner0 + icoord;
      int col1 = 3*indCorner3 + icoord;
      int col2 = 3*indVertex0 + icoord;
      int col3 = 3*indVertex1 + icoord;

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
    nrow += 3;
  }
  nvertex = (int)sideVertexIndices[5].size()-ny;
  for (ivertex=1; ivertex<nvertex; ivertex++)
  {
     if (ivertex%ny==ny-1) 
      continue;

    int indVertex0 =  sideVertexIndices[4][ivertex];
    int indVertex1 =  sideVertexIndices[5][ivertex];

    for (int icoord=0; icoord<3; icoord++)
    {
      int row = nrow+icoord;

      int col0 = 3*indCorner0 + icoord;
      int col1 = 3*indCorner5 + icoord;
      int col2 = 3*indVertex0 + icoord;
      int col3 = 3*indVertex1 + icoord;

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
    nrow += 3;
  }

  for(int ii = 0;ii<nconstraints;ii++){
    triplets.push_back(TripletS(nrowInit+ii, nrowInit+ii,0.f));
  }
  KConstraint.setFromTriplets(triplets.begin(), triplets.end());
  K += KConstraint; 
}




