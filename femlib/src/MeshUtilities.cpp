#include "MeshUtilities.h"

#include "ElementRegGrid.hpp"
#include "Element.hpp"
#include "ElementRegGrid2D.h"
#include "Element2D.h"
#include "Material2D.h"
#include "Material.hpp"
#include "cfgMaterialUtilities.h"
#include <Eigen/Dense>
#include "Tensor.h"
#include "ElementHex2D.h"
#include "ElementHex.hpp"
#include "StrainLin2D.h"
#include "StrainLin.hpp"
#include "MaterialQuad2D.h"
#include "Quadrature2D.h"
#include "cfgUtilities.h"

int meshUtil::getElementIndex(int nx, int ny, int nz, int ii , int jj, int kk)
{
  return ii * ny * nz + jj * nz + kk;
}

//0: Left, 1: Right, 2: Bottom, 3:Top, 4: Back, 5: Front
void meshUtil::getSideVoxels(int iSide, int nx, int ny, int nz, std::vector<int> &oSideVoxelIndices)
{
  oSideVoxelIndices.clear();

   int ii, jj, kk;
   if (iSide==0) //left
   {    
     for(jj=0; jj<ny; jj++)
     {
       for(kk=0; kk<nz; kk++)
       {
         int indElement = meshUtil::getElementIndex(nx, ny, nz, 0, jj, kk);
         oSideVoxelIndices.push_back(indElement);
       }
     }
   }
   else if (iSide==1) //right
   {
     for(jj=0; jj<ny; jj++)
     {
       for(kk=0; kk<nz; kk++)
       {
         int indElement = meshUtil::getElementIndex(nx, ny, nz, nx-1, jj, kk);
         oSideVoxelIndices.push_back(indElement);
       }
     }
   }
   else if (iSide==2) //bottom
   {
     for(ii=0; ii<nx; ii++)
     {
       for(kk=0; kk<nz; kk++)
       {
         int indElement = meshUtil::getElementIndex(nx, ny, nz, ii, 0, kk);
         oSideVoxelIndices.push_back(indElement);
       }
     }
   }
   else if (iSide==3) //top
   {
     for(ii=0; ii<nx; ii++)
     {
       for(kk=0; kk<nz; kk++)
       {
         int indElement = meshUtil::getElementIndex(nx, ny, nz, ii, ny-1, kk);
         oSideVoxelIndices.push_back(indElement);
       }
     }
   }
   else if (iSide==4) //back
   {
     for(ii=0; ii<nx; ii++)
     {
       for(jj=0; jj<ny; jj++)
       {
         int indElement = meshUtil::getElementIndex(nx, ny, nz, ii, jj, 0);
         oSideVoxelIndices.push_back(indElement);
       }
     } 
   }
   else if (iSide==5) //front
   {
     for(ii=0; ii<nx; ii++)
     {
       for(jj=0; jj<ny; jj++)
       {
         int indElement = meshUtil::getElementIndex(nx, ny, nz, ii, jj, nz-1);
         oSideVoxelIndices.push_back(indElement);
       }
     }
   }
}

//0: Left, 1: Right, 2: Bottom, 3:Top, 4: Back, 5: Front
void meshUtil::getSideVertices(int iSide, const ElementRegGrid * iElementGrid, std::vector<int> &oElemIndices, std::vector<std::vector<int> > &oFaceVertexIndices)
{
  assert(iElementGrid);
  
  oElemIndices.clear();
  oFaceVertexIndices.clear();

  int nx = iElementGrid->nx;
  int ny = iElementGrid->ny;
  int nz = iElementGrid->nz;

   // [axis][side][face_vertex]
   int externalFaceIndices[3][2][4] = {
    { {0,1,2,3}, {4,5,6,7} },
    { {0,1,4,5}, {2,3,6,7} }, 
    { {0,4,2,6}, {1,5,3,7} } };

   int ii, jj, kk;
   if (iSide==0) //left
   {    
     for(jj=0; jj<ny; jj++)
     {
       for(kk=0; kk<nz; kk++)
       {
         int indElement = iElementGrid->GetEleInd(0,jj,kk);
         oElemIndices.push_back(indElement);
         std::vector<int> fvIndices;
         for (int ivertex=0; ivertex<4; ivertex++)
         {
           fvIndices.push_back(externalFaceIndices[0][0][ivertex]);
         }
         oFaceVertexIndices.push_back(fvIndices);
       }
     }
   }
   else if (iSide==1) //right
   {
     for(jj=0; jj<ny; jj++)
     {
       for(kk=0; kk<nz; kk++)
       {
         int indElement = iElementGrid->GetEleInd(nx-1,jj,kk);
         oElemIndices.push_back(indElement);
         std::vector<int> fvIndices;
         for (int ivertex=0; ivertex<4; ivertex++)
         {
           fvIndices.push_back(externalFaceIndices[0][1][ivertex]);
         }
         oFaceVertexIndices.push_back(fvIndices);
       }
     }
   }
   else if (iSide==2) //bottom
   {
     for(ii=0; ii<nx; ii++)
     {
       for(kk=0; kk<nz; kk++)
       {
         int indElement = iElementGrid->GetEleInd(ii,0,kk);
         oElemIndices.push_back(indElement);
         std::vector<int> fvIndices;
         for (int ivertex=0; ivertex<4; ivertex++)
         {
           fvIndices.push_back(externalFaceIndices[1][0][ivertex]);
         }
         oFaceVertexIndices.push_back(fvIndices);
       }
     }
   }
   else if (iSide==3) //top
   {
     for(ii=0; ii<nx; ii++)
     {
       for(kk=0; kk<nz; kk++)
       {
         int indElement = iElementGrid->GetEleInd(ii,ny-1,kk);
         oElemIndices.push_back(indElement);
         std::vector<int> fvIndices;
         for (int ivertex=0; ivertex<4; ivertex++)
         {
           fvIndices.push_back(externalFaceIndices[1][1][ivertex]);
         }
         oFaceVertexIndices.push_back(fvIndices);
       }
     }
   }
   else if (iSide==4) //back
   {
     for(ii=0; ii<nx; ii++)
     {
       for(jj=0; jj<ny; jj++)
       {
         int indElement = iElementGrid->GetEleInd(ii,jj,0);
         oElemIndices.push_back(indElement);
         std::vector<int> fvIndices;
         for (int ivertex=0; ivertex<4; ivertex++)
         {
           fvIndices.push_back(externalFaceIndices[2][0][ivertex]);
         }
         oFaceVertexIndices.push_back(fvIndices);
       }
     } 
   }
   else if (iSide==5) //front
   {
     for(ii=0; ii<nx; ii++)
     {
       for(jj=0; jj<ny; jj++)
       {
         int indElement = iElementGrid->GetEleInd(ii,jj,nz-1);
         oElemIndices.push_back(indElement);
         std::vector<int> fvIndices;
         for (int ivertex=0; ivertex<4; ivertex++)
         {
           fvIndices.push_back(externalFaceIndices[2][1][ivertex]);
         }
         oFaceVertexIndices.push_back(fvIndices);
       }
     }
   }
}

//0: Left, 1: Right, 2: Bottom, 3:Top, 4: Back, 5: Front
void meshUtil::getSideVertices(int iSide, const ElementRegGrid * iElementGrid, std::vector<int> &oVertexIndices)
{
  assert(iElementGrid);
  
  oVertexIndices.clear();

  int nx = iElementGrid->nx;
  int ny = iElementGrid->ny;
  int nz = iElementGrid->nz;

   // [axis][side][face_vertex]
   int externalFaceIndices[3][2][4] = {
    { {0,1,2,3}, {4,5,6,7} },
    { {0,4,1,5}, {2,6,3,7} }, 
    { {0,2,4,6}, {1,3,5,7} } };

   int ii, jj, kk;
   if (iSide==0) //left
   {    
     for(jj=0; jj<ny; jj++)
     {
       for(kk=0; kk<nz; kk++)
       {
         int indElement = iElementGrid->GetEleInd(0,jj,kk);

         int fvIndex = externalFaceIndices[0][0][0];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
         if (kk==nz-1)
         {
           int fvIndex = externalFaceIndices[0][0][1];
           int indVertex = iElementGrid->e[indElement]->at(fvIndex);
           oVertexIndices.push_back(indVertex);
         }
       }
     }
     for(kk=0; kk<nz; kk++)
     {
       int indElement = iElementGrid->GetEleInd(0,ny-1,kk);

       int fvIndex = externalFaceIndices[0][0][2];
       int indVertex = iElementGrid->e[indElement]->at(fvIndex);
       oVertexIndices.push_back(indVertex);
       if (kk==nz-1)
       {
         int fvIndex = externalFaceIndices[0][0][3];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
       }
     }
   }
   else if (iSide==1) //right
   {
     for(jj=0; jj<ny; jj++)
     {
       for(kk=0; kk<nz; kk++)
       {
         int indElement = iElementGrid->GetEleInd(nx-1,jj,kk);

         int fvIndex = externalFaceIndices[0][1][0];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
         if (kk==nz-1)
         {
           int fvIndex = externalFaceIndices[0][1][1];
           int indVertex = iElementGrid->e[indElement]->at(fvIndex);
           oVertexIndices.push_back(indVertex);
         }
       }
     }
     for(kk=0; kk<nz; kk++)
     {
       int indElement = iElementGrid->GetEleInd(nx-1,ny-1,kk);

       int fvIndex = externalFaceIndices[0][1][2];
       int indVertex = iElementGrid->e[indElement]->at(fvIndex);
       oVertexIndices.push_back(indVertex);
       if (kk==nz-1)
       {
         int fvIndex = externalFaceIndices[0][1][3];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
       }
     }
   }
   else if (iSide==2) //bottom
   {
     for(kk=0; kk<nz; kk++)
     {
       for(ii=0; ii<nx; ii++)
       {
         int indElement = iElementGrid->GetEleInd(ii,0,kk);

         int fvIndex = externalFaceIndices[1][0][0];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
          if (ii==nx-1)
         {
           int fvIndex = externalFaceIndices[1][0][1];
           int indVertex = iElementGrid->e[indElement]->at(fvIndex);
           oVertexIndices.push_back(indVertex);
         }
       }
     }
     for(ii=0; ii<nx; ii++)
     {
       int indElement = iElementGrid->GetEleInd(ii,0,nz-1);

       int fvIndex = externalFaceIndices[1][0][2];
       int indVertex = iElementGrid->e[indElement]->at(fvIndex);
       oVertexIndices.push_back(indVertex);
       if (ii==nx-1)
       {
         int fvIndex = externalFaceIndices[1][0][3];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
       }
     }
   }
   else if (iSide==3) //top
   {
     for(kk=0; kk<nz; kk++) 
     {
       for(ii=0; ii<nx; ii++)
       {
         int indElement = iElementGrid->GetEleInd(ii,ny-1,kk);

         int fvIndex = externalFaceIndices[1][1][0];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
         if (ii==nx-1)
         {
           int fvIndex = externalFaceIndices[1][1][1];
           int indVertex = iElementGrid->e[indElement]->at(fvIndex);
           oVertexIndices.push_back(indVertex);
         }
       }
     }
     for(ii=0; ii<nx; ii++)
     {
       int indElement = iElementGrid->GetEleInd(ii,ny-1,nz-1);

       int fvIndex = externalFaceIndices[1][1][2];
       int indVertex = iElementGrid->e[indElement]->at(fvIndex);
       oVertexIndices.push_back(indVertex);
       if (ii==nx-1)
       {
         int fvIndex = externalFaceIndices[1][1][3];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
       }
     }
   }
   else if (iSide==4) //back
   {
     for(ii=0; ii<nx; ii++)
     {
       for(jj=0; jj<ny; jj++)
       {
         int indElement = iElementGrid->GetEleInd(ii,jj,0);

         int fvIndex = externalFaceIndices[2][0][0];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
         if (jj==ny-1)
         {
           int fvIndex = externalFaceIndices[2][0][1];
           int indVertex = iElementGrid->e[indElement]->at(fvIndex);
           oVertexIndices.push_back(indVertex);
         }
       }
     } 
     for(jj=0; jj<ny; jj++)
     {
       int indElement = iElementGrid->GetEleInd(nx-1,jj,0);

       int fvIndex = externalFaceIndices[2][0][2];
       int indVertex = iElementGrid->e[indElement]->at(fvIndex);
       oVertexIndices.push_back(indVertex);
       if (jj==ny-1)
       {
         int fvIndex = externalFaceIndices[2][0][3];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
       }
     }
   }
   else if (iSide==5) //front
   {
     for(ii=0; ii<nx; ii++)
     {
       for(jj=0; jj<ny; jj++)
       {
         int indElement = iElementGrid->GetEleInd(ii,jj,nz-1);

         int fvIndex = externalFaceIndices[2][1][0];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
         if (jj==ny-1)
         {
           int fvIndex = externalFaceIndices[2][1][1];
           int indVertex = iElementGrid->e[indElement]->at(fvIndex);
           oVertexIndices.push_back(indVertex);
         }
       }
     }
     for(jj=0; jj<ny; jj++)
     {
       int indElement = iElementGrid->GetEleInd(nx-1,jj,nz-1);

       int fvIndex = externalFaceIndices[2][1][2];
       int indVertex = iElementGrid->e[indElement]->at(fvIndex);
       oVertexIndices.push_back(indVertex);
       if (jj==ny-1)
       {
         int fvIndex = externalFaceIndices[2][1][3];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
       }
     }
   }
}

//0: Left, 1: Right, 2: Bottom, 3:Top
void meshUtil::getSideVertices(int iSide, const ElementRegGrid2D * iElementGrid, std::vector<int> &oElemIndices, std::vector<std::vector<int> > &oEdgeVertexIndices)
{
  assert(iElementGrid);
  
  oElemIndices.clear();
  oEdgeVertexIndices.clear();

  int nx = iElementGrid->nx;
  int ny = iElementGrid->ny;

  // [axis][side][edge_vertex]
   int externalFaceIndices[2][2][2] = {
    { {0,1}, {2,3} },
    { {0,2}, {1,3} } };

   int ii, jj;
   if (iSide==0) //left
   {    
     for(jj=0; jj<ny; jj++)
     {
         int indElement = iElementGrid->GetEleInd(0,jj);
         oElemIndices.push_back(indElement);
         std::vector<int> evIndices;
         for (int ivertex=0; ivertex<2; ivertex++)
         {
           evIndices.push_back(externalFaceIndices[0][0][ivertex]);
         }
         oEdgeVertexIndices.push_back(evIndices);
     }
   }
   else if (iSide==1) //right
   {
     for(jj=0; jj<ny; jj++)
     {
       int indElement = iElementGrid->GetEleInd(nx-1,jj);
       oElemIndices.push_back(indElement);
       std::vector<int> evIndices;
       for (int ivertex=0; ivertex<2; ivertex++)
       {
         evIndices.push_back(externalFaceIndices[0][1][ivertex]);
       }
       oEdgeVertexIndices.push_back(evIndices);
     }
   }
   else if (iSide==2) //bottom
   {
     for(ii=0; ii<nx; ii++)
     {
       int indElement = iElementGrid->GetEleInd(ii,0);
       oElemIndices.push_back(indElement);
       std::vector<int> evIndices;
       for (int ivertex=0; ivertex<2; ivertex++)
       {
         evIndices.push_back(externalFaceIndices[1][0][ivertex]);
       }
       oEdgeVertexIndices.push_back(evIndices);
     }
   }
   else if (iSide==3) //top
   {
     for(ii=0; ii<nx; ii++)
     {
       int indElement = iElementGrid->GetEleInd(ii,ny-1);
       oElemIndices.push_back(indElement);
       std::vector<int> evIndices;
       for (int ivertex=0; ivertex<2; ivertex++)
       {
         evIndices.push_back(externalFaceIndices[1][1][ivertex]);
       }
       oEdgeVertexIndices.push_back(evIndices);
     }
   }
}

//0: Left, 1: Right, 2: Bottom, 3:Top
void meshUtil::getSideVertices(int iSide, const ElementRegGrid2D * iElementGrid, std::vector<int> &oVertexIndices)
{
  assert(iElementGrid);
  
  oVertexIndices.clear();

  int nx = iElementGrid->nx;
  int ny = iElementGrid->ny;

  // [axis][side][edge_vertex]
   int externalFaceIndices[2][2][2] = {
    { {0,1}, {2,3} },
    { {0,2}, {1,3} } };

   int ii, jj;
   if (iSide==0) //left
   {    
     for(jj=0; jj<ny; jj++)
     {
       int indElement = iElementGrid->GetEleInd(0,jj);
       if (jj==0)
       {
         int fvIndex = externalFaceIndices[0][0][0];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
       }
       int fvIndex = externalFaceIndices[0][0][1];
       int indVertex = iElementGrid->e[indElement]->at(fvIndex);
       oVertexIndices.push_back(indVertex);
     }
   }
   else if (iSide==1) //right
   {
     for(jj=0; jj<ny; jj++)
     {
       int indElement = iElementGrid->GetEleInd(nx-1,jj);
       if (jj==0)
       {
         int fvIndex = externalFaceIndices[0][1][0];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
       }
       int fvIndex = externalFaceIndices[0][1][1];
       int indVertex = iElementGrid->e[indElement]->at(fvIndex);
       oVertexIndices.push_back(indVertex);
     }
   }
   else if (iSide==2) //bottom
   {
     for(ii=0; ii<nx; ii++)
     {
       int indElement = iElementGrid->GetEleInd(ii,0);
       if (ii==0)
       {
         int fvIndex = externalFaceIndices[1][0][0];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
       }
       int fvIndex = externalFaceIndices[1][0][1];
       int indVertex = iElementGrid->e[indElement]->at(fvIndex);
       oVertexIndices.push_back(indVertex);
     }
   }
   else if (iSide==3) //top
   {
     for(ii=0; ii<nx; ii++)
     {
       int indElement = iElementGrid->GetEleInd(ii,ny-1);
       if (ii==0)
       {
         int fvIndex = externalFaceIndices[1][1][0];
         int indVertex = iElementGrid->e[indElement]->at(fvIndex);
         oVertexIndices.push_back(indVertex);
       }
       int fvIndex = externalFaceIndices[1][1][1];
       int indVertex = iElementGrid->e[indElement]->at(fvIndex);
       oVertexIndices.push_back(indVertex);
     }
   }
}

void meshUtil::getExternalVertices(ElementRegGrid * iElementGrid, std::vector<int> &oElemIndices, std::vector<std::vector<int> > &oFaceVertexIndices, std::vector<Vector3S> &oNormals)
{
  assert(iElementGrid);

  oElemIndices.clear();
  oFaceVertexIndices.clear();
  oNormals.clear();

  int nx = iElementGrid->nx;
  int ny = iElementGrid->ny;
  int nz = iElementGrid->nz;

  int externalFaceIndices[3][2][4] = {
    { {0,1,2,3}, {4,5,6,7} },
    { {0,1,4,5}, {2,3,6,7} }, 
    { {0,4,2,6}, {1,5,3,7} } };

  int ii, jj, kk;
  //left
  for(jj=0; jj<ny; jj++)
  {
    for(kk=0; kk<nz; kk++)
    {
      int indElement = iElementGrid->GetEleInd(0,jj,kk);
      oElemIndices.push_back(indElement);
      std::vector<int> fvIndices;
      for (int ivertex=0; ivertex<4; ivertex++)
      {
        fvIndices.push_back(externalFaceIndices[0][0][ivertex]);
      }
      oFaceVertexIndices.push_back(fvIndices);
      oNormals.push_back(Vector3S(-1, 0, 0));
    }
  }
  //right
  for(jj=0; jj<ny; jj++)
  {
    for(kk=0; kk<nz; kk++)
    {
      int indElement = iElementGrid->GetEleInd(nx-1,jj,kk);
      oElemIndices.push_back(indElement);
      std::vector<int> fvIndices;
      for (int ivertex=0; ivertex<4; ivertex++)
      {
        fvIndices.push_back(externalFaceIndices[0][1][ivertex]);
      }
      oFaceVertexIndices.push_back(fvIndices);
      oNormals.push_back(Vector3S(1, 0, 0));
    }
  }
  //bottom
  for(ii=0; ii<nx; ii++)
  {
    for(kk=0; kk<nz; kk++)
    {
      int indElement = iElementGrid->GetEleInd(ii,0,kk);
      oElemIndices.push_back(indElement);
      std::vector<int> fvIndices;
      for (int ivertex=0; ivertex<4; ivertex++)
      {
        fvIndices.push_back(externalFaceIndices[1][0][ivertex]);
      }
      oFaceVertexIndices.push_back(fvIndices);
      oNormals.push_back(Vector3S(0, -1, 0));
    }
  }
  //top
  for(ii=0; ii<nx; ii++)
  {
    for(kk=0; kk<nz; kk++)
    {
      int indElement = iElementGrid->GetEleInd(ii,ny-1,kk);
      oElemIndices.push_back(indElement);
      std::vector<int> fvIndices;
      for (int ivertex=0; ivertex<4; ivertex++)
      {
        fvIndices.push_back(externalFaceIndices[1][1][ivertex]);
      }
      oFaceVertexIndices.push_back(fvIndices);
      oNormals.push_back(Vector3S(0, 1, 0));
    }
  }
  //back
  for(ii=0; ii<nx; ii++)
  {
    for(jj=0; jj<ny; jj++)
    {
      int indElement = iElementGrid->GetEleInd(ii,jj,0);
      oElemIndices.push_back(indElement);
      std::vector<int> fvIndices;
      for (int ivertex=0; ivertex<4; ivertex++)
      {
        fvIndices.push_back(externalFaceIndices[2][0][ivertex]);
      }
      oFaceVertexIndices.push_back(fvIndices);
      oNormals.push_back(Vector3S(0, 0, -1));
    }
  } 
  //front
  for(ii=0; ii<nx; ii++)
  {
    for(jj=0; jj<ny; jj++)
    {
      int indElement = iElementGrid->GetEleInd(ii,jj,nz-1);
      oElemIndices.push_back(indElement);
      std::vector<int> fvIndices;
      for (int ivertex=0; ivertex<4; ivertex++)
      {
        fvIndices.push_back(externalFaceIndices[2][1][ivertex]);
      }
      oFaceVertexIndices.push_back(fvIndices);
      oNormals.push_back(Vector3S(0, 0, 1));
    }
  }
}

void meshUtil::getVertexValences(const ElementRegGrid * iElementGrid, const std::vector<int> &iElemIndices, const std::vector<std::vector<int> > &iFvIndices, std::map<int,int> &oind2Valence)
{
  oind2Valence.clear();
  int ielem=0, nelem=(int)iElemIndices.size();
  for (ielem=0; ielem<nelem; ielem++)
  {
    int indElement = iElemIndices[ielem];
    int ivertex=0, nvertex=(int)iFvIndices[ielem].size();
    for (ivertex=0; ivertex<nvertex; ivertex++)
    {
      int fvIndex = iFvIndices[ielem][ivertex];
      int indVertex = iElementGrid->e[indElement]->at(fvIndex);
      oind2Valence[indVertex]++;
    }
  }
}

void meshUtil::getVertexValences(const ElementRegGrid2D * iElementGrid, const std::vector<int> &iElemIndices, const std::vector<std::vector<int> > &iEvIndices, std::map<int,int> &oind2Valence)
{
  oind2Valence.clear();
  int ielem=0, nelem=(int)iElemIndices.size();
  for (ielem=0; ielem<nelem; ielem++)
  {
    int indElement = iElemIndices[ielem];
    int ivertex=0, nvertex=(int)iEvIndices[ielem].size();
    for (ivertex=0; ivertex<nvertex; ivertex++)
    {
      int fvIndex = iEvIndices[ielem][ivertex];
      int indVertex = iElementGrid->e[indElement]->at(fvIndex);
      oind2Valence[indVertex]++;
    }
  }
}

cfgScalar meshUtil::getVonMisesStress(int iElementIndex, ElementMesh2D &iElementGrid)
{
  cfgScalar maxStress = 0;

  int matIndex = iElementGrid.me[iElementIndex];
  Material2D * mat = iElementGrid.m[matIndex];
  assert(mat);
  Element2D * element = iElementGrid.e[iElementIndex];
  assert(element);
  std::vector<Matrix2S> stressTensors = mat->getStressTensors(element, &iElementGrid);
  int itensor=0, ntensor=(int)stressTensors.size();
  for (itensor=0; itensor<ntensor; itensor++)
  {
    Matrix2S & tensor = stressTensors[itensor];

    std::vector<cfgScalar> stresses;
    stresses.push_back(tensor(0,0));
    stresses.push_back(tensor(1,1));
    stresses.push_back(tensor(1,0));

    cfgScalar stress = cfgMaterialUtilities::computeVonMisesStress2D(stresses);
    if (stress > maxStress)
    {
      maxStress = stress;
    }
  }
  return maxStress;
 }

std::vector<cfgScalar> meshUtil::getVonMisesStressesPerElement(ElementMesh2D &iElementGrid)
{
  std::vector<cfgScalar> stresses;
  int ielem, nelem=(int)iElementGrid.e.size();
  for (ielem=0; ielem<nelem; ielem++)
  {
    cfgScalar stress = getVonMisesStress(ielem, iElementGrid);
    stresses.push_back(stress);
  }
  return stresses;
 }

void meshUtil::computeStrains(ElementRegGrid2D &iElementGrid, const std::vector<std::vector<cfgScalar> > &iHarmonicDisplacements, std::vector<std::vector<cfgScalar> > &oStrains)
{
  oStrains.clear();
  oStrains.resize(iHarmonicDisplacements.size());

  int nx = iElementGrid.nx;
  int ny = iElementGrid.ny;

  std::vector<int> corners;
  corners.push_back(iElementGrid.GetVertInd(0, 0));
  corners.push_back(iElementGrid.GetVertInd(0, ny));
  corners.push_back(iElementGrid.GetVertInd(nx, 0));
  corners.push_back(iElementGrid.GetVertInd(nx, ny));
  ElementHex2D coarseElem(corners);
  StrainLin2D strainLin;
  
  int idisp, ndisp=(int)iHarmonicDisplacements.size();
  for (idisp=0; idisp<ndisp; idisp++)
  {
    std::vector<Vector2S> h = cfgMaterialUtilities::toVector2S(iHarmonicDisplacements[idisp]);
    std::vector<Vector2S> x = cfgUtil::add(iElementGrid.X, h);
    Matrix2S coarseF = coarseElem.defGrad(Vector2S(0,0), iElementGrid.X, x);
    Matrix2S currentStrainTensor = strainLin.getStrainTensor(coarseF);
    oStrains[idisp].push_back(currentStrainTensor(0,0));
    oStrains[idisp].push_back(currentStrainTensor(1,1));
    oStrains[idisp].push_back(currentStrainTensor(0,1));
  }
}

void meshUtil::computeStrains(ElementRegGrid &iElementGrid, const std::vector<std::vector<cfgScalar> > &iHarmonicDisplacements, std::vector<std::vector<cfgScalar> > &oStrains)
{
  oStrains.clear();
  oStrains.resize(iHarmonicDisplacements.size());

  int nx = iElementGrid.nx;
  int ny = iElementGrid.ny;
  int nz = iElementGrid.nz;

  std::vector<int> corners;
  corners.push_back(iElementGrid.GetVertInd(0, 0, 0));
  corners.push_back(iElementGrid.GetVertInd(0, 0, nz));
  corners.push_back(iElementGrid.GetVertInd(0, ny, 0));
  corners.push_back(iElementGrid.GetVertInd(0, ny, nz));
  corners.push_back(iElementGrid.GetVertInd(nx, 0, 0));
  corners.push_back(iElementGrid.GetVertInd(nx, 0, nz));
  corners.push_back(iElementGrid.GetVertInd(nx, ny, 0));
  corners.push_back(iElementGrid.GetVertInd(nx, ny, nz));
  ElementHex coarseElem(corners);
  StrainLin strainLin;

  int idisp, ndisp=(int)iHarmonicDisplacements.size();
  assert(ndisp==6);
  for (idisp=0; idisp<ndisp; idisp++)
  {
    std::vector<Eigen::Vector3f> h = cfgMaterialUtilities::toVector3f(cfgUtil::convertVec<cfgScalar,float>(iHarmonicDisplacements[idisp]));
    std::vector<Eigen::Vector3f> x = cfgUtil::add(iElementGrid.X, h);
    Eigen::Matrix3f coarseF = coarseElem.defGrad(Eigen::Vector3f::Zero(), iElementGrid.X, x);
    Eigen::Matrix3f currentStrainTensor = strainLin.getStrainTensor(coarseF);
    oStrains[idisp].push_back(currentStrainTensor(0,0));
    oStrains[idisp].push_back(currentStrainTensor(1,1));
    oStrains[idisp].push_back(currentStrainTensor(2,2));
    oStrains[idisp].push_back(currentStrainTensor(1,2));
    oStrains[idisp].push_back(currentStrainTensor(0,2));
    oStrains[idisp].push_back(currentStrainTensor(0,1));
  }
}

void meshUtil::computeCoarsenedElasticityTensor(ElementRegGrid2D &iElementGrid, const std::vector<std::vector<cfgScalar> > &iHarmonicDisplacements, MatrixXS &oCoarsenedTensor)
{
  int nx = iElementGrid.nx;
  int ny = iElementGrid.ny;

  std::vector<int> corners;
  corners.push_back(iElementGrid.GetVertInd(0, 0));
  corners.push_back(iElementGrid.GetVertInd(0, ny));
  corners.push_back(iElementGrid.GetVertInd(nx, 0));
  corners.push_back(iElementGrid.GetVertInd(nx, ny));
  ElementHex2D coarseElem(corners);
  StrainLin2D strainLin;
  
  std::vector<Matrix2S> coarseStrainTensors;
  int idisp, ndisp=(int)iHarmonicDisplacements.size();
  assert(ndisp==3);
  for (idisp=0; idisp<ndisp; idisp++)
  {
    std::vector<Vector2S> h = cfgMaterialUtilities::toVector2S(iHarmonicDisplacements[idisp]);
    std::vector<Vector2S> x = cfgUtil::add(iElementGrid.X, h);
    Matrix2S coarseF = coarseElem.defGrad(Vector2S(0,0), iElementGrid.X, x);
    Matrix2S currentStrainTensor = strainLin.getStrainTensor(coarseF);
    coarseStrainTensors.push_back(currentStrainTensor);
    //std::cout << "strain = " << currentStrainTensor(0,0) << " " << currentStrainTensor(0,1) << " " << currentStrainTensor(1,1) << std::endl;
  }

  /*std::vector<std::vector<Matrix2S> > coarseG(2);
  coarseG[0].push_back(coarseStrainTensors[0]);
  coarseG[0].push_back(coarseStrainTensors[2]);
  coarseG[1].push_back(coarseStrainTensors[2]);
  coarseG[1].push_back(coarseStrainTensors[1]);

  MatrixXS coarseGv = Tensor::toVoigtRepresentation(coarseG);*/ 

  MatrixXS coarseGv = MatrixXS::Zero(3,3);
  for (idisp=0; idisp<ndisp; idisp++)
  {
    coarseGv(0, idisp) = coarseStrainTensors[idisp](0,0);
    coarseGv(1, idisp) = coarseStrainTensors[idisp](1,1);
    coarseGv(2, idisp) = coarseStrainTensors[idisp](0,1);
  }

  MatrixXS invCoarseGv = coarseGv.inverse();
  //std::cout << "G = " << coarseGv << std::endl;

  MatrixXS Gt_C_G = MatrixXS::Zero(3, 3);

  int ielem, nelem=(int)iElementGrid.e.size();
  for (ielem=0; ielem<nelem; ielem++)
  {
    int matIndex = iElementGrid.me[ielem];
    Material2D * mat = iElementGrid.m[matIndex];
    assert(mat);

    Element2D * element = iElementGrid.e[ielem];
    assert(element);

    std::vector<MatrixXS> C = mat->getElasticityTensors();

    std::vector<std::vector<Matrix2S> > strainTensors;
    int idisp, ndisp=(int)iHarmonicDisplacements.size();
    assert(ndisp==3);
    for (idisp=0; idisp<ndisp; idisp++)
    {
      std::vector<Vector2S> h = cfgMaterialUtilities::toVector2S(iHarmonicDisplacements[idisp]);
      std::vector<Vector2S> x = cfgUtil::add(iElementGrid.X, h);
      std::vector<Matrix2S> currentStrainTensors = mat->getStrainTensors(element, &iElementGrid, x);
      strainTensors.push_back(currentStrainTensors);
    }

    std::vector<cfgScalar> w = ((MaterialQuad2D*)mat)->q->w;

    int itensor, ntensor=(int)strainTensors[0].size();
    for (itensor=0; itensor<ntensor; itensor++)
    {
      /*std::vector<std::vector<Matrix2S> > G(2);
      G[0].push_back(strainTensors[0][itensor]);
      G[0].push_back(strainTensors[2][itensor]);
      G[1].push_back(strainTensors[2][itensor]);
      G[1].push_back(strainTensors[1][itensor]);

      MatrixXS Gv = Tensor::toVoigtRepresentation(G);*/ 

      MatrixXS Gv = MatrixXS::Zero(3,3);
      for (idisp=0; idisp<ndisp; idisp++)
      {
        Gv(0, idisp) = strainTensors[idisp][itensor](0,0);
        Gv(1, idisp) = strainTensors[idisp][itensor](1,1);
        Gv(2, idisp) = strainTensors[idisp][itensor](0,1);
      }

      Gt_C_G += w[itensor]*Gv.transpose()*C[itensor]*Gv;

      //std::cout << "Gv = " << Gv << std::endl;
    }
  }
  oCoarsenedTensor = (1.f/nelem) * invCoarseGv.transpose() * Gt_C_G * invCoarseGv;
}

void meshUtil::computeCoarsenedElasticityTensor(ElementRegGrid &iElementGrid, const std::vector<std::vector<cfgScalar> > &iHarmonicDisplacements, MatrixXS &oCoarsenedTensor)
{
  int nx = iElementGrid.nx;
  int ny = iElementGrid.ny;
  int nz = iElementGrid.nz;

  std::vector<int> corners;
  corners.push_back(iElementGrid.GetVertInd(0, 0, 0));
  corners.push_back(iElementGrid.GetVertInd(0, 0, nz));
  corners.push_back(iElementGrid.GetVertInd(0, ny, 0));
  corners.push_back(iElementGrid.GetVertInd(0, ny, nz));
  corners.push_back(iElementGrid.GetVertInd(nx, 0, 0));
  corners.push_back(iElementGrid.GetVertInd(nx, 0, nz));
  corners.push_back(iElementGrid.GetVertInd(nx, ny, 0));
  corners.push_back(iElementGrid.GetVertInd(nx, ny, nz));
  ElementHex coarseElem(corners);
  StrainLin strainLin;
  
  std::vector<Eigen::Matrix3f> coarseStrainTensors;
  int idisp, ndisp=(int)iHarmonicDisplacements.size();
  assert(ndisp==6);
  for (idisp=0; idisp<ndisp; idisp++)
  {
    std::vector<Eigen::Vector3f> h = cfgMaterialUtilities::toVector3f(cfgUtil::convertVec<cfgScalar, float>(iHarmonicDisplacements[idisp]));
    std::vector<Eigen::Vector3f> x = cfgUtil::add(iElementGrid.X, h);
    Eigen::Matrix3f coarseF = coarseElem.defGrad(Eigen::Vector3f::Zero(), iElementGrid.X, x);
    Eigen::Matrix3f currentStrainTensor = strainLin.getStrainTensor(coarseF);
    coarseStrainTensors.push_back(currentStrainTensor);
    //std::cout << "strain = " << currentStrainTensor(0,0) << " " << currentStrainTensor(0,1) << " " << currentStrainTensor(1,1) << std::endl;
  }

  MatrixXS coarseGv = MatrixXS::Zero(6,6);
  for (idisp=0; idisp<ndisp; idisp++)
  {
    coarseGv(0, idisp) = coarseStrainTensors[idisp](0,0);
    coarseGv(1, idisp) = coarseStrainTensors[idisp](1,1);
    coarseGv(2, idisp) = coarseStrainTensors[idisp](2,2);
    coarseGv(3, idisp) = coarseStrainTensors[idisp](1,2);
    coarseGv(4, idisp) = coarseStrainTensors[idisp](0,2);
    coarseGv(5, idisp) = coarseStrainTensors[idisp](0,1);
  }

  MatrixXS invCoarseGv = coarseGv.inverse();
  //std::cout << "G = " << coarseGv << std::endl;

  MatrixXS Gt_C_G = MatrixXS::Zero(6, 6);

  int ielem, nelem=(int)iElementGrid.e.size();
  for (ielem=0; ielem<nelem; ielem++)
  {
    int matIndex = iElementGrid.me[ielem];
    Material * mat = iElementGrid.m[matIndex];
    assert(mat);

    Element* element = iElementGrid.e[ielem];
    assert(element);

    std::vector<Eigen::MatrixXf> Ctmp = mat->getElasticityTensors();
    std::vector<MatrixXS> C;
    for (int imat=0; imat<Ctmp.size(); imat++)
    {
      C.push_back(cfgMaterialUtilities::toMatrixScalar(Ctmp[imat]));
    }
    if (0)
    {
      std::cout << "C = " <<  C[0] << std::endl;
    }

    std::vector<std::vector<Eigen::Matrix3f> > strainTensors;
    int idisp, ndisp=(int)iHarmonicDisplacements.size();
    assert(ndisp==6);
    for (idisp=0; idisp<ndisp; idisp++)
    {
      std::vector<Eigen::Vector3f> h = cfgMaterialUtilities::toVector3f(cfgUtil::convertVec<cfgScalar, float>(iHarmonicDisplacements[idisp]));
      std::vector<Eigen::Vector3f> x = cfgUtil::add(iElementGrid.X, h);
      std::vector<Eigen::Matrix3f> currentStrainTensors = mat->getStrainTensors(element, &iElementGrid, x);
      strainTensors.push_back(currentStrainTensors);
    }

    std::vector<cfgScalar> w = ((MaterialQuad2D*)mat)->q->w;

    int itensor, ntensor=(int)strainTensors[0].size();
    for (itensor=0; itensor<ntensor; itensor++)
    {
      MatrixXS Gv = MatrixXS::Zero(6,6);
      for (idisp=0; idisp<ndisp; idisp++)
      {
        Gv(0, idisp) = strainTensors[idisp][itensor](0,0);
        Gv(1, idisp) = strainTensors[idisp][itensor](1,1);
        Gv(2, idisp) = strainTensors[idisp][itensor](2,2);
        Gv(3, idisp) = strainTensors[idisp][itensor](1,2);
        Gv(4, idisp) = strainTensors[idisp][itensor](0,2);
        Gv(5, idisp) = strainTensors[idisp][itensor](0,1); 
      }
      Gt_C_G += w[itensor]*Gv.transpose()*C[itensor]*Gv;

      //std::cout << "Gv = " << Gv << std::endl;
    }
  }
  oCoarsenedTensor = (1.f/nelem) * invCoarseGv.transpose() * Gt_C_G * invCoarseGv;
}






