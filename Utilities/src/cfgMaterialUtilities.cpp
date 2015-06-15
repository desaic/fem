#include "cfgMaterialUtilities.h"

#include <iostream>
#include <fstream>
#include <assert.h>

#include "ElementRegGrid.hpp"
#include "Element.hpp"

#include "qhullUtilities.h"

#include <Qhull.h>
#include <QhullFacetList.h>
#include <QhullVertex.h>
#include <QhullVertexSet.h>
using namespace orgQhull;

#include <cfgUtilities.h>
using namespace cfgUtil;

bool cfgMaterialUtilities::readMaterialCombinations(const std::string iFileName, std::vector<std::vector<int> > &oMaterials)
{
  std::ifstream stream(iFileName);
  if (!stream.is_open())
  {
    return false;
  }

  deSerialize<int>(stream, oMaterials, "Materials");
  return true;
}


bool cfgMaterialUtilities::writeMaterialCombinations(const std::string iFileName, const std::vector<std::vector<int> > &iMaterials)
{
  std::ofstream stream(iFileName);
  if (!stream.is_open())
  {
    std::cout << "cfgMaterialUtilities::writeMaterialCombinations: invalid fileName" << std::endl;
    return false;
  }
   
  stream << "Materials" << std::endl;
  int ivec=0, nvec=(int)iMaterials.size();
  stream << nvec << " " << std::endl;
  for (ivec=0; ivec<nvec; ivec++)
  {
    const std::vector<int> &iMaterialVec = iMaterials[ivec];
    int imat=0, nmat=(int)iMaterialVec.size();
    stream << nmat << " ";
    for (imat=0; imat<nmat; imat++)
    {
      stream << iMaterialVec[imat] << " ";
    }
    stream << std::endl;
  }
  return true;
}

bool cfgMaterialUtilities::readData(const std::string &iFileName,  Vector3f &oForceAxis, std::vector<int> &oMaterials, std::vector<float> &oStresses, std::vector<float> &oStrains)
{
  std::ifstream stream(iFileName);
  if (!stream.is_open())
  {
    return false;
  }

  std::string dummy;
  // Force Axis
  stream >> dummy; 
  stream >> oForceAxis[0] >>  oForceAxis[1] >> oForceAxis[2];
  
  // Materials
  deSerialize<int>(stream, oMaterials,"Materials");

  // Stress strain
  std::vector<float> StrainStress[2];
  deSerialize<float,2>(stream, StrainStress, "StressStrain");
  oStrains = StrainStress[0];
  oStresses = StrainStress[1];

  return true;
}

bool cfgMaterialUtilities::writeData(const std::string &iFileName, const Vector3f &iForceAxis, const std::vector<int> &iMaterials, const std::vector<float> &iStresses, const std::vector<float> &iStrains)
{
  std::ofstream stream(iFileName);
  if (!stream.is_open())
  {
    return false;
  }

  stream << "ForceAxis" << std::endl; 
  stream << iForceAxis[0] << " " << iForceAxis[1] << " " << iForceAxis[2]  << std::endl;
  
  stream << "Materials" << std::endl;
  int imat=0, nmat=(int)iMaterials.size();
  stream << nmat << " ";
  for (imat=0; imat<nmat; imat++)
  {
    stream << iMaterials[imat] << " ";
  }
  stream << std::endl;

  stream << "StressStrain" << std::endl;
  int isample, nsample=(int)iStresses.size();
  stream << nsample+1 << " " << std::endl;
  stream << 0 << " " << 0 << std::endl; 
  for (isample=0; isample<nsample; isample++)
  {
    stream << iStrains[isample] << " " << iStresses[isample]  << std::endl; 
  }
  return false;
}

cfgScalar cfgMaterialUtilities::computeYoungModulus(const std::vector<cfgScalar> &iStrains, const std::vector<cfgScalar>  &iStresses)
{
  cfgScalar YoungModulus=0;
  fitLinearFunction<cfgScalar>(iStrains, iStresses, YoungModulus);
  return YoungModulus;
}

bool cfgMaterialUtilities::computeMaterialParameters(const std::string &iMaterialFile, const std::string iStressStrainFilesDirectories[2], const std::string iStressStrainBaseFileName, int iNbFiles,
                                                    std::vector<cfgScalar> &oPhysicalParameters, std::vector<std::vector<int> > &oMaterialAssignments)
{
  std::vector<std::vector<int> > baseMaterialStructures;
  readMaterialCombinations(iMaterialFile, baseMaterialStructures);

  oPhysicalParameters.clear(); //YoungModulus, density;
  int icomb, ncomb=iNbFiles;
  oMaterialAssignments.clear();
  oMaterialAssignments.resize(ncomb);
  int naxis = 2;
  for (icomb=0; icomb<ncomb; icomb++)
  {
    int iaxis;
    for (iaxis=0; iaxis<naxis; iaxis++)
    {
      std::string sampleFile = iStressStrainFilesDirectories[iaxis] + iStressStrainBaseFileName + "_" + std::to_string(icomb) + ".txt"; 
      Vector3f forceAxis;
      std::vector<int> materials; 
      std::vector<cfgScalar> stresses, strains;
      readData(sampleFile, forceAxis, materials, stresses, strains);

      if (iaxis==0)
      {
        cfgScalar r=0;
        int imat, nmat=(int)materials.size();
        for (imat=0; imat<nmat; imat++)
        {
          int matIndex = materials[imat];
          std::vector<int> & baseMaterialStructure = baseMaterialStructures[matIndex];
          int icell, ncell=(int)baseMaterialStructure.size();
          for (icell=0; icell<ncell; icell++)
          {
            r += baseMaterialStructure[icell];
          }
          oMaterialAssignments[icomb].insert(oMaterialAssignments[icomb].end(), baseMaterialStructure.begin(), baseMaterialStructure.end());
        }
        int ncell = (int)baseMaterialStructures[0].size()*nmat;
        r /= (cfgScalar)ncell;
        oPhysicalParameters.push_back(r);
      }
      cfgScalar Y = computeYoungModulus(strains, stresses);
      oPhysicalParameters.push_back(Y);
    }
  }
  return true;
}

void cfgMaterialUtilities::getMaterialAssignment(int nx, int ny, int nz, const std::vector<int> &iMaterialCombIndices, int nX, int nY, int nZ, const std::vector<std::vector<int> > &iBaseCellMaterials, std::vector<int> &oMaterials)
{
  oMaterials.clear();
  oMaterials.resize(nx*ny*nz, -1);

  int repX = nx/(2*nX);
  int repY = ny/(2*nY);
  int repZ = nz/(2*nZ);

  int i,j,k;
  for (i=0; i<repX; i++)
  {
    for (j=0; j<repY; j++)
    {
      for (k=0; k<repZ; k++)
      {
        int l, m, n;
        for (l=0; l<2; l++)
        {
          for (m=0; m<2; m++)
          {
            for (n=0; n<2; n++)
            {
              int matCombinationIndex = iMaterialCombIndices[l* 4 + m * 2 + n];
              const std::vector<int> &matCombination = iBaseCellMaterials[matCombinationIndex];

              int ii,jj,kk;
              for (ii=0; ii<nX; ii++)
              {
                for (jj=0; jj<nY; jj++)
                {
                  for (kk=0; kk<nZ; kk++)
                  {
                    int indx = i*(2*nX) + nX*l + ii;
                    int indy = j*(2*nY) + nY*m + jj;
                    int indz = k*(2*nZ) + nZ*n + kk;
                    int elementIndex = indx * ny*nz + indy * nz + indz;

                    oMaterials[elementIndex] = matCombination[ii*nY*nZ + jj*nZ + kk];
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

int cfgMaterialUtilities::getElementIndex(int nx, int ny, int nz, int ii , int jj, int kk)
{
  return ii * ny * nz + jj * nz + kk;
}

//0: Left, 1: Right, 2: Bottom, 3:Top, 4: Back, 5: Front
void cfgMaterialUtilities::getSideVoxels(int iSide, int nx, int ny, int nz, std::vector<int> &oSideVoxelIndices)
{
  oSideVoxelIndices.clear();

   int ii, jj, kk;
   if (iSide==0) //left
   {    
     for(jj=0; jj<ny; jj++)
     {
       for(kk=0; kk<nz; kk++)
       {
         int indElement = cfgMaterialUtilities::getElementIndex(nx, ny, nz, 0, jj, kk);
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
         int indElement = cfgMaterialUtilities::getElementIndex(nx, ny, nz, nx-1, jj, kk);
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
         int indElement = cfgMaterialUtilities::getElementIndex(nx, ny, nz, ii, 0, kk);
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
         int indElement = cfgMaterialUtilities::getElementIndex(nx, ny, nz, ii, ny-1, kk);
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
         int indElement = cfgMaterialUtilities::getElementIndex(nx, ny, nz, ii, jj, 0);
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
         int indElement = cfgMaterialUtilities::getElementIndex(nx, ny, nz, ii, jj, nz-1);
         oSideVoxelIndices.push_back(indElement);
       }
     }
   }
}

//0: Left, 1: Right, 2: Bottom, 3:Top, 4: Back, 5: Front
void cfgMaterialUtilities::getSideVertices(int iSide, const ElementRegGrid * iElementGrid, std::vector<int> &oElemIndices, std::vector<std::vector<int> > &oFaceVertexIndices)
{
  assert(iElementGrid);
  
  oElemIndices.clear();
  oFaceVertexIndices.clear();

  int nx = iElementGrid->nx;
  int ny = iElementGrid->ny;
  int nz = iElementGrid->nz;

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

void cfgMaterialUtilities::getExternalVertices(ElementRegGrid * iElementGrid, std::vector<int> &oElemIndices, std::vector<std::vector<int> > &oFaceVertexIndices, std::vector<Vector3f> &oNormals)
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
      oNormals.push_back(Vector3f(-1.f, 0.f, 0.f));
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
      oNormals.push_back(Vector3f(1.f, 0.f, 0.f));
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
      oNormals.push_back(Vector3f(0.f, -1.f, 0.f));
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
      oNormals.push_back(Vector3f(0.f, 1.f, 0.f));
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
      oNormals.push_back(Vector3f(0.f, 0.f, -1.f));
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
      oNormals.push_back(Vector3f(0.f, 0.f, 1.f));
    }
  }
}

void cfgMaterialUtilities::getVertexValences(const ElementRegGrid * iElementGrid, const std::vector<int> &iElemIndices, const std::vector<std::vector<int> > &iFvIndices, std::map<int,int> &oind2Valence)
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


void cfgMaterialUtilities::computeConvexHull(const std::vector<float> &iPoints, int iDim, std::vector<int> &oConvexHullVertices)
{
  assert(iPoints.size()%iDim==0);

  oConvexHullVertices.clear();
  coordT * points = cfgUtil::createArray<float,coordT>(iPoints);
  int numpoints = (int)iPoints.size()/iDim;
  boolT newIsMalloc = 0;

  /*qh_init_A (stdin, stdout, stderr, 0, NULL);
  int exitcode = setjmp (qh errexit);
  if (!exitcode)
  {
    qh_init_B(points, numpoints, iDim, newIsMalloc);
    qh_qhull();
    qh_check_output();
    print_summary();
  }*/ 

  std::string rboxCommand2, qhullCommand2;
  Qhull qhull(rboxCommand2.c_str(), iDim, numpoints, points, qhullCommand2.c_str());

  QhullVertexList vertices = qhull.vertexList();
  QhullLinkedList<QhullVertex>::const_iterator v_it, v_end=vertices.end();
  for (v_it=vertices.begin(); v_it!=v_end; ++v_it)
  {
    QhullPoint p = (*v_it).point();
    oConvexHullVertices.push_back(p.id());

    /*const coordT *coord = p.coordinates();
    int icoord;
    for (icoord=0; icoord<iDim; icoord++)
    {
      oConvexHullVertices.push_back(coord[icoord]);
    }*/ 
  }

  /*qhull.outputQhull();
  if(qhull.hasQhullMessage()){
    std::cout << "\nResults of qhull\n" << qhull.qhullMessage();
    qhull.clearQhullMessage();
  }

  QhullFacetList facets= qhull.facetList();
  std::cout << "\nFacets created by Qhull::runQhull()\n" << facets;

  QhullVertexList vertices = qhull.vertexList();
  std::cout << "\nVertices created by Qhull::runQhull()\n" << vertices; */ 
}

void cfgMaterialUtilities::computeDelaundayTriangulation(const std::vector<float> &iPoints, int iDim, std::vector<int> &oFaceVertices, std::vector<int> &oBoundaryVertices, std::vector<int> &oBoundaryFaces)
{
  assert(iPoints.size()%iDim==0);

  oFaceVertices.clear();
  coordT * points = cfgUtil::createArray<float,coordT>(iPoints);
  int numpoints = (int)iPoints.size()/iDim;
  boolT newIsMalloc = 0;

  std::string rboxCommand2, qhullCommand2="d QJ";
  Qhull qhull(rboxCommand2.c_str(), iDim, numpoints, points, qhullCommand2.c_str());

  std::vector<int> allFaceVertices;
  qhullUtilities::getFacetsVertices(qhull.facetList(), allFaceVertices);
  //oFaceVertices = allFaceVertices;


  float meanEdgeLength = getMeanEdgeLength(iPoints, allFaceVertices, 4);
  float minEdgeLength = getMinEdgeLength(iPoints, allFaceVertices, 4);
  //qhullUtilities::getBoundaryVertices(qhull, qhull.facetList(), 1.65, *oBoundaryVertices);

  std::vector<QhullFacet> smallFacets, largeFacets;
  qhullUtilities::sortFaces(qhull, qhull.facetList(), 1.8, smallFacets, largeFacets);
  qhullUtilities::getBoundaryVertices(smallFacets, oBoundaryVertices, oBoundaryFaces);

  oFaceVertices.clear();
  qhullUtilities::getFacetsVertices(smallFacets, oFaceVertices);

  /*oFaceVertices.clear();
  qhullUtilities::filterFaces(qhull, qhull.facetList(), 1.5, oFaceVertices);
  *oBoundaryVertices = oFaceVertices;*/ 


  if (0)
  {
    std::ofstream stream("Summary.txt");
    qhull.outputQhull();
    if(qhull.hasQhullMessage()){
      stream << "\nResults of qhull\n" << qhull.qhullMessage();
      qhull.clearQhullMessage();
    }

    QhullFacetList facets= qhull.facetList();
    stream << "\nFacets created by Qhull::runQhull()\n" << facets;

    QhullVertexList vertices = qhull.vertexList();
    stream << "\nVertices created by Qhull::runQhull()\n" << vertices;
  }
}

float cfgMaterialUtilities::getMeanEdgeLength(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim)
{
  assert(iIndexArray.size()%iDim==0);

  float MeanEdgeLength=0;
  int nedge=0;
  int iface=0, nface=(int)iIndexArray.size()/iDim;
  for (iface=0; iface<nface; iface++)
  {
    int ivertex=0;
    for (ivertex=0; ivertex<iDim; ivertex++)
    {
      int indVertex1 = iIndexArray[iDim*iface+ivertex];
      int indVertex2 = iIndexArray[iDim*iface+(ivertex+1)%iDim];

      Vector3f p1 = getVector3f(indVertex1, iX);
      Vector3f p2 = getVector3f(indVertex2, iX);
      float length = (p2-p1).abs();
      MeanEdgeLength += length;
      nedge++;
    }
  }
  if (nedge>0)
    MeanEdgeLength /= nedge;

  return MeanEdgeLength;
}

float cfgMaterialUtilities::getMinEdgeLength(const std::vector<float> &iX,  const std::vector<int> &iIndexArray, int iDim)
{
  assert(iIndexArray.size()%iDim==0);

  float SqMinEdgeLength=FLT_MAX;
  int iface=0, nface=(int)iIndexArray.size()/iDim;
  for (iface=0; iface<nface; iface++)
  {
    int ivertex=0;
    for (ivertex=0; ivertex<iDim; ivertex++)
    {
      int indVertex1 = iIndexArray[iDim*iface+ivertex];
      int indVertex2 = iIndexArray[iDim*iface+(ivertex+1)%iDim];

      Vector3f p1 = getVector3f(indVertex1, iX);
      Vector3f p2 = getVector3f(indVertex2, iX);
      float SqLength = (p2-p1).absSquared();
      if (SqLength < SqMinEdgeLength)
      {
        SqMinEdgeLength = SqLength;
      }
    }
  }
  return sqrt(SqMinEdgeLength);
}

void cfgMaterialUtilities::getEdgesFromTriFaceIndexArray(const std::vector<int> &iTriIndexArray, std::vector<int> &oEdgeIndexArray)
{
  std::vector<int> edgeIndices;
  int ifacet=0, nfacet=(int)iTriIndexArray.size()/3;
  for (ifacet=0; ifacet<nfacet; ifacet++)
  {
    int ipoint=0;
    for (ipoint=0; ipoint<3; ipoint++)
    {
      int indPoint1 = iTriIndexArray[3*ifacet+ipoint];
      oEdgeIndexArray.push_back(indPoint1);

      int indPoint2 = iTriIndexArray[3*ifacet+(ipoint+1)%3];
      oEdgeIndexArray.push_back(indPoint2);
    }
  }
}

void cfgMaterialUtilities::getEdgesFromTetFaceIndexArray(const std::vector<int> &iTetIndexArray, std::vector<int> &oEdgeIndexArray)
{
  oEdgeIndexArray.clear();

  int tetfaces[4][3] = { {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 3, 0} };

  std::vector<int> edgeIndices;
  int ifacet=0, nfacet=(int)iTetIndexArray.size()/4;
  for (ifacet=0; ifacet<nfacet; ifacet++)
  {
    int iface=0;
    for (iface=0; iface<4; iface++)
    {
      int ipoint=0;
      for (ipoint=0; ipoint<4; ipoint++)
      {
        int indVertex1 = tetfaces[iface][ipoint%3];
        int indPoint1 = iTetIndexArray[4*ifacet+indVertex1];
        oEdgeIndexArray.push_back(indPoint1);

        int indVertex2 = tetfaces[iface][(ipoint+1)%3];
        int indPoint2 = iTetIndexArray[4*ifacet+indVertex2];
        oEdgeIndexArray.push_back(indPoint2);
      }
    }
  }
}

Vector3f cfgMaterialUtilities::getVector3f(int indVertex, const std::vector<float> &iPoints)
{
  assert(iPoints.size()%3==0);
  int indPoint = 3*indVertex;
  return Vector3f(iPoints[indPoint], iPoints[indPoint+1], iPoints[indPoint+2]);
}
