
#include "SamplesGeneratorImpl.h"

#include <assert.h>

#include "ElementRegGrid.hpp"
#include "ElementRegGrid2D.h"
#include "Element.hpp"
#include "Element2D.h"

#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"
#include "StrainLin.hpp"

#include "MaterialQuad2D.h"
#include "StrainLin2D.h"

#include "Stepper.hpp"
#include "StepperNewton.hpp"
#include "Stepper2D.h"
#include "StepperNewton2D.h"
#include "AdmmCPU.hpp"
#include "StepperGrad.hpp"

#include <algorithm>

#include "cfgUtilities.h"
using namespace cfgUtil;

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;


SamplesGeneratorImpl::SamplesGeneratorImpl(int iDim)
{
  m_UseLinearMaterial = true;
  assert(iDim==2||iDim==3);
  m_dim = iDim;
  m_blockRep = 3;
}

SamplesGeneratorImpl::~SamplesGeneratorImpl()
{
}

float SamplesGeneratorImpl::computeStrain(const ElementRegGrid * iElementGrid, int iAxis)
{
  //writeVector2File(iElementGrid->x, m_OutputDirectory + "x.m");

  Vector3f Barycenters[2];

  int iside;
  for (iside=0; iside<2; iside++)
  {
    std::vector<int> elemIndices;
    std::vector<std::vector<int> > fvIndices;
    getSideVertices(2*iAxis+iside, iElementGrid, elemIndices, fvIndices);

    std::map<int,int> ind2Valence;
    getVertexValences(iElementGrid, elemIndices, fvIndices, ind2Valence);

    
    int ielem=0, nelem=(int)elemIndices.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      int indElement = elemIndices[ielem];
      int ivertex=0, nvertex=(int)fvIndices[ielem].size();
      for (ivertex=0; ivertex<nvertex; ivertex++)
      {
        int fvIndex = fvIndices[ielem][ivertex];
        int indVertex = iElementGrid->e[indElement]->at(fvIndex);
        float weight = 1.f/ind2Valence[indVertex];
        Barycenters[iside] += weight*iElementGrid->x[indVertex];
      }
    }
    Barycenters[iside] /= (float)ind2Valence.size();
  }
  float strain = (Barycenters[0]-Barycenters[1]).abs()-1;
  return strain;
}

float SamplesGeneratorImpl::computeStrain(const ElementRegGrid2D * iElementGrid, int iAxis)
{
  //writeVector2File(iElementGrid->x, m_OutputDirectory + "x.m");

  Vector2f Barycenters[2];

  int iside;
  for (iside=0; iside<2; iside++)
  {
    std::vector<int> elemIndices;
    std::vector<std::vector<int> > evIndices;
    getSideVertices(2*iAxis+iside, iElementGrid, elemIndices, evIndices);

    std::map<int,int> ind2Valence;
    getVertexValences(iElementGrid, elemIndices, evIndices, ind2Valence);

    
    int ielem=0, nelem=(int)elemIndices.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      int indElement = elemIndices[ielem];
      int ivertex=0, nvertex=(int)evIndices[ielem].size();
      for (ivertex=0; ivertex<nvertex; ivertex++)
      {
        int fvIndex = evIndices[ielem][ivertex];
        int indVertex = iElementGrid->e[indElement]->at(fvIndex);
        float weight = 1.f/ind2Valence[indVertex];
        Barycenters[iside] += weight*iElementGrid->x[indVertex];
      }
    }
    Barycenters[iside] /= (float)ind2Valence.size();
  }
  float strain = (Barycenters[0]-Barycenters[1]).abs()-1;
  return strain;
}


void SamplesGeneratorImpl::setExternalForces(ElementRegGrid * iElementGrid, int iAxis, int iSide, float iForceMagnitude)
{
  assert(iElementGrid);

  ElementRegGrid * em = iElementGrid;
  std::vector<int> elemIndices;
  std::vector<std::vector<int> > fvIndices;
  getSideVertices(2*iAxis+iSide, em, elemIndices, fvIndices);

  int n[3];
  n[0] = em->nx;
  n[1] = em->ny;
  n[2] = em->nz;

  Vector3f ff;
  ff[iAxis] = (1.f/(n[(iAxis+1)%3]*n[(iAxis+2)%3])) * iForceMagnitude;

  int ielem=0, nelem=(int)elemIndices.size();
  for (ielem=0; ielem<nelem; ielem++)
  {
    int indElement = elemIndices[ielem];
    int ivertex=0, nvertex=(int)fvIndices[ielem].size();
    for (ivertex=0; ivertex<nvertex; ivertex++)
    {
      int FvIndex = fvIndices[ielem][ivertex];
      int VertexIndex =em->e[indElement]->at(FvIndex);
      em->fe[VertexIndex] += ff;
    }
  }

  // Constraints
  std::vector<int> elemIndicesOppositeSide;
  std::vector<std::vector<int> > fvIndicesOppositeSide;
  getSideVertices(2*iAxis+1-iSide, em, elemIndicesOppositeSide, fvIndicesOppositeSide);
  int nelemOppositeSide=(int)elemIndicesOppositeSide.size();
  for (ielem=0; ielem<nelemOppositeSide; ielem++)
  {
    int indElement = elemIndicesOppositeSide[ielem];
    int ivertex=0, nvertex=(int)fvIndicesOppositeSide[ielem].size();
    for (ivertex=0; ivertex<nvertex; ivertex++)
    {
      int FvIndex = fvIndicesOppositeSide[ielem][ivertex];
      int VertexIndex =em->e[indElement]->at(FvIndex);
      em->fixed[VertexIndex] = 1;
    }
  }
}

void SamplesGeneratorImpl::setExternalForces(ElementRegGrid2D * iElementGrid, int iAxis, int iSide, float iForceMagnitude)
{
  assert(iElementGrid);

  ElementRegGrid2D * em = iElementGrid;
  std::vector<int> elemIndices;
  std::vector<std::vector<int> > fvIndices;
  getSideVertices(2*iAxis+iSide, em, elemIndices, fvIndices);

  int n[2];
  n[0] = em->nx;
  n[1] = em->ny;

  double thickness = 1.;

  Vector2f ff;
  ff[iAxis] = (1.f/(n[(iAxis+1)%2]*thickness)) * iForceMagnitude;

  int ielem=0, nelem=(int)elemIndices.size();
  for (ielem=0; ielem<nelem; ielem++)
  {
    int indElement = elemIndices[ielem];
    int ivertex=0, nvertex=(int)fvIndices[ielem].size();
    for (ivertex=0; ivertex<nvertex; ivertex++)
    {
      int FvIndex = fvIndices[ielem][ivertex];
      int VertexIndex =em->e[indElement]->at(FvIndex);
      em->fe[VertexIndex] += ff;
    }
  }

  // Constraints
  std::vector<int> elemIndicesOppositeSide;
  std::vector<std::vector<int> > fvIndicesOppositeSide;
  getSideVertices(2*iAxis+1-iSide, em, elemIndicesOppositeSide, fvIndicesOppositeSide);
  int nelemOppositeSide=(int)elemIndicesOppositeSide.size();
  for (ielem=0; ielem<nelemOppositeSide; ielem++)
  {
    int indElement = elemIndicesOppositeSide[ielem];
    int ivertex=0, nvertex=(int)fvIndicesOppositeSide[ielem].size();
    for (ivertex=0; ivertex<nvertex; ivertex++)
    {
      int FvIndex = fvIndicesOppositeSide[ielem][ivertex];
      int VertexIndex =em->e[indElement]->at(FvIndex);
      em->fixed[VertexIndex] = 1;
    }
  }
}

ElementRegGrid * SamplesGeneratorImpl::createPhysicalSystem(int iN[3], std::vector<MaterialQuad> &iMaterials)
{
  ElementRegGrid * em = new ElementRegGrid(iN[0],iN[1],iN[2]);

  // Materials
  int imat=0, nmat=(int)iMaterials.size();
  for (imat=0; imat<nmat; imat++)
  {
    em->addMaterial(&iMaterials[imat]);
  }
  em->check();
  return em;
}

ElementRegGrid2D * SamplesGeneratorImpl::createPhysicalSystem(int iN[2], std::vector<MaterialQuad2D> &iMaterials)
{
  ElementRegGrid2D * em = new ElementRegGrid2D(iN[0],iN[1]);

  // Materials
  int imat=0, nmat=(int)iMaterials.size();
  for (imat=0; imat<nmat; imat++)
  {
    em->addMaterial(&iMaterials[imat]);
  }
  em->check();
  return em;
}

/*
ElementRegGrid * SamplesGeneratorImpl::createPhysicalSystem(int iN[3], Vector3f iForce, std::vector<MaterialQuad> &iMaterials)
{
  ElementRegGrid * em = new ElementRegGrid(iN[0],iN[1],iN[2]);

  std::vector<int> elemIndices;
  std::vector<std::vector<int> > fvIndices;
  std::vector<Vector3f> Normals;
  getExternalVertices(em, elemIndices, fvIndices, Normals);

  // External forces
  Vector3f ff;
  int iaxis;
  for (iaxis=0; iaxis<3; iaxis++)
  {
    ff[iaxis] = (1.f/(iN[(iaxis+1)%3]*iN[(iaxis+2)%3])) * iForce[iaxis];
  }

  int ielem=0, nelem=(int)fvIndices.size();
  for (ielem=0; ielem<nelem; ielem++)
  {
    int indElement = elemIndices[ielem];
    Vector3f & Normal = Normals[ielem];
    std::vector<int> & elemFvIndices = fvIndices[ielem];
    if (Vector3f::dot(ff,Normal)>0)
    {
      Vector3f Force = fabs(Vector3f::dot(ff,Normal))*Normal;
      int ivertex=0, nvertex=(int)elemFvIndices.size();
      for (ivertex=0; ivertex<nvertex; ivertex++)
      {
        int FvIndex = elemFvIndices[ivertex];
        int VertexIndex =em->e[indElement]->at(FvIndex);
        em->fe[VertexIndex] += Force;
      }
    }
  }

  // Constraints
   for(int ii = 0;ii<iN[0];ii++){
    for(int jj =0;jj<iN[2];jj++){
      int eidx= em->GetEleInd(ii,0,jj);
      int aa[4] = {0,1,4,5};
      for(int kk = 0;kk<4;kk++){
        int vidx =em->e[eidx]->at(aa[kk]);
        em->fixed[vidx] = 1;
      }
    }
  }

  // Materials
  int imat=0, nmat=(int)iMaterials.size();
  for (imat=0; imat<nmat; imat++)
  {
    em->addMaterial(&iMaterials[imat]);
  }
  em->check();
  return em;
}*/ 

bool SamplesGeneratorImpl::computeEquilibrium(Stepper * iStepper, ElementRegGrid * iElementGrid, int iNumberOfSteps, bool &oConverged)
{
  assert(iStepper);

  iStepper->nSteps = iNumberOfSteps;
  iStepper->init(iElementGrid);
  int status = iStepper->step();
  oConverged = status<0;
  return status<=0;
}

bool SamplesGeneratorImpl::computeEquilibrium(Stepper2D * iStepper, ElementRegGrid2D * iElementGrid, int iNumberOfSteps, bool &oConverged)
{
  assert(iStepper);

  iStepper->nSteps = iNumberOfSteps;
  iStepper->init(iElementGrid);
  int status = iStepper->step();
  oConverged = status<0;
  return status<=0;
}

bool SamplesGeneratorImpl::computeForceDeformationSample(ElementRegGrid * iElementGrid, std::string iStepperType, int iNumberOfSteps, bool &oConverged)
{
  Stepper * stepper = 0;
  if (iStepperType == "newton"){
    stepper = new StepperNewton();
  }
  else if (iStepperType == "newtonCuda"){
 //   stepper = new NewtonCuda();
 //   ((NewtonCuda*)stepper)->m_fixRigidMotion = false;
  }
  else if (iStepperType == "ipopt"){
    //not yet implemented
    //  stepper = new IpoptStepper();
  }
  else if (iStepperType == "admmCPU"){
    stepper = new AdmmCPU();
    ((AdmmCPU*)stepper)->ro0 = 10;
  }
  else if (iStepperType == "grad"){
    stepper = new StepperGrad();
  }
  else{
    //stepper = new AdmmNoSpring();
    //stepper = new ADMMStepper();  
  }
  //stepper->rmRigid = true;
  bool Converged = false;
  bool ResOk = computeEquilibrium(stepper, iElementGrid, iNumberOfSteps, oConverged);
  
  delete(stepper);

  return ResOk;
}

bool SamplesGeneratorImpl::computeForceDeformationSample(ElementRegGrid2D * iElementGrid, std::string iStepperType, int iNumberOfSteps, bool &oConverged)
{
  Stepper2D * stepper = 0;
  if (iStepperType == "newton"){
    stepper = new StepperNewton2D();
  }
  else if (iStepperType == "newtonCuda"){
 //   stepper = new NewtonCuda();
 //   ((NewtonCuda*)stepper)->m_fixRigidMotion = false;
  }
  else if (iStepperType == "ipopt"){
    //not yet implemented
    //  stepper = new IpoptStepper();
  }
  else if (iStepperType == "admmCPU"){
   // stepper = new AdmmCPU();
   // ((AdmmCPU*)stepper)->ro0 = 10;
  }
  else if (iStepperType == "grad"){
  //  stepper = new StepperGrad();
  }
  else{
    //stepper = new AdmmNoSpring();
    //stepper = new ADMMStepper();  
  }
  //stepper->rmRigid = true;
  bool Converged = false;
  bool ResOk = computeEquilibrium(stepper, iElementGrid, iNumberOfSteps, oConverged);
  
  delete(stepper);

  return ResOk;
}

void SamplesGeneratorImpl::assignMaterials(ElementRegGrid * iElementGrid, int iMacroStructureMaterials[2][2][2])
{
  assert(iElementGrid);

  int nx = iElementGrid->nx;
  int ny = iElementGrid->ny;
  int nz = iElementGrid->nz;

  int i,j,k;
  for (i=0; i<nx; i++)
  {
    for (j=0; j<ny; j++)
    {
      for (k=0; k<nz; k++)
      {
        int materialIndex = iMacroStructureMaterials[i%2][j%2][k%2];
        int elementIndex = iElementGrid->GetEleInd(i, j, k);
        iElementGrid->me[elementIndex] = materialIndex;
      }
    }
  }
}

int SamplesGeneratorImpl::checkBorderSimilarity(const std::vector<int> &iMaterialAssignment, int N[3], const std::vector<std::vector<int> > &iBaseMaterialStructures)
{
  int NbMismatches = 0;

  //0: Left, 1: Right, 2: Bottom, 3:Top, 4: Back, 5: Front

  int externalFaceIndices[3][2][4] = {
    { {0,1,2,3}, {4,5,6,7} },
    { {0,1,4,5}, {2,3,6,7} }, 
    { {0,4,2,6}, {1,5,3,7} } };

   int iaxis;
   for (iaxis=0; iaxis<3; iaxis++)
   {
     std::vector<int> voxels1, voxels2;
     cfgMaterialUtilities::getSideVoxels(2*iaxis, N[0], N[1], N[2], voxels1);
     cfgMaterialUtilities::getSideVoxels(2*iaxis+1, N[0], N[1], N[2], voxels2);

     int iface;
     for (iface=0; iface<4; iface++)
     {
       int iside, nside=2;
       for (iside=0; iside<nside; iside++)
       {
         int indFace1 = externalFaceIndices[iaxis][iside][iface];
         int indFace2 = externalFaceIndices[iaxis][1-iside][iface];

         int indMatStructure1 = iMaterialAssignment[indFace1];
         int indMatStructure2 = iMaterialAssignment[indFace2];
         if (indMatStructure1>=0 && indMatStructure2>=0)
         {
           std::vector<int> materials1 = iBaseMaterialStructures[indMatStructure1];
           std::vector<int> materials2 = iBaseMaterialStructures[indMatStructure2];
           int ivox, nvox=(int)voxels1.size();
           for (ivox=0; ivox<nvox; ivox++)
           {
             int indVoxel1 = voxels1[ivox];
             int indMat1 =  materials1[indVoxel1];

             int indVoxel2 = voxels2[ivox];
             int indMat2 =  materials2[indVoxel2];

             if (indMat1 != indMat2)
             {
               NbMismatches++;
             }
           }
         }
       }
     }
   }
  return NbMismatches;
}

void SamplesGeneratorImpl::sampleMaterialSpace()
{
  int level = 2;
  int blockSize = level+1;

  std::string materialFile = m_OutputDirectory + "Materials_" + std::to_string(level) + ".txt";
  std::vector<std::vector<int> > materialCombinations;
  readMaterialCombinations(materialFile, materialCombinations);

  std::vector<cfgScalar> physicalParameters; //YoungModulus, density;
  std::vector<std::vector<int> > macroCellMaterials;
  int icomb, ncomb=256;
  macroCellMaterials.resize(ncomb);
  for (icomb=0; icomb<ncomb; icomb++)
  {
    std::string sampleFile = m_OutputDirectory + "StressStrain_" + std::to_string(level) + "_" + std::to_string(blockSize) + "_" + std::to_string(icomb) + ".txt"; 
    Vector3f forceAxis;
    std::vector<int> materials;
    std::vector<cfgScalar> stresses, strains;
    readData(sampleFile, forceAxis, materials, stresses, strains);

    cfgScalar r=0;
    int imat, nmat=(int)materials.size();
    for (imat=0; imat<nmat; imat++)
    {
      int matIndex = materials[imat];
      std::vector<int> & voxelMaterials = materialCombinations[matIndex];
      int icell, ncell=(int)voxelMaterials.size();
      for (icell=0; icell<ncell; icell++)
      {
        r += voxelMaterials[icell];
      }
      macroCellMaterials[icomb].insert(macroCellMaterials[icomb].end(), voxelMaterials.begin(), voxelMaterials.end());
    }
    int ncell = materialCombinations[0].size()*nmat;
    r /= (cfgScalar)ncell;

    cfgScalar Y = computeYoungModulus(strains, stresses);
    physicalParameters.push_back(Y);
    physicalParameters.push_back(r);
  }
  std::vector<int> convexHullVertices;
  computeConvexHull(physicalParameters, 2, convexHullVertices);
  std::vector<std::vector<int> > newMaterialCombinations = getSubVector(macroCellMaterials, convexHullVertices);
}

bool SamplesGeneratorImpl::sampleDeformation(int iN[3], std::vector<MaterialQuad> &iMaterial, const std::string iStepperType, float iMaxForce, int iForceAxis, int iNumberOfSample, const std::vector<int> & iMaterials,
                                             std::vector<float> &oStresses, std::vector<float> &oStrains)
{
  bool ResOk = true;

  oStresses.clear();
  oStrains.clear();

  std::vector<Vector3f> xinit;
  int isample;
  for (isample=0; isample<iNumberOfSample && ResOk; isample++)
  {
    float forceMagnitude = (isample+1.f)*iMaxForce/iNumberOfSample;

    ElementRegGrid * em = createPhysicalSystem(iN, iMaterial);
    setExternalForces(em, iForceAxis, 1, forceMagnitude);
    em->me = iMaterials;
    //assignMaterials(em, iMaterials);

    if (isample>0)
    {
      em->x = xinit;
    }
    bool Converged = false;
    ResOk = computeForceDeformationSample(em, iStepperType, 100, Converged);
    if (ResOk)
    {
      float strain = computeStrain(em, iForceAxis);
      xinit = em->x;

      //writeVector2File(em->me, m_OutputDirectory + "materials.m");
      std::cout << "Strain = " << strain << std::endl;

      delete em;

      oStresses.push_back(forceMagnitude);
      oStrains.push_back(strain);
    }
  }
  return ResOk;
}

bool SamplesGeneratorImpl::sampleDeformation(int iN[2], std::vector<MaterialQuad2D> &iMaterial, const std::string iStepperType, float iMaxForce, int iForceAxis, int iNumberOfSample, const std::vector<int> & iMaterials,
                                             std::vector<float> &oStresses, std::vector<float> &oStrains)
{
  bool ResOk = true;

  oStresses.clear();
  oStrains.clear();

  std::vector<Vector2f> xinit;
  int isample;
  for (isample=0; isample<iNumberOfSample && ResOk; isample++)
  {
    float forceMagnitude = (isample+1.f)*iMaxForce/iNumberOfSample;

    ElementRegGrid2D * em = createPhysicalSystem(iN, iMaterial);
    setExternalForces(em, iForceAxis, 1, forceMagnitude);
    em->me = iMaterials;
    //assignMaterials(em, iMaterials);

    if (isample>0)
    {
      em->x = xinit;
    }
    bool Converged = false;
    ResOk = computeForceDeformationSample(em, iStepperType, 100, Converged);
    if (ResOk)
    {
      float strain = computeStrain(em, iForceAxis);
      xinit = em->x;

      //writeVector2File(em->me, m_OutputDirectory + "materials.m");
      std::cout << "Strain = " << strain << std::endl;

      delete em;

      oStresses.push_back(forceMagnitude);
      oStrains.push_back(strain);
    }
  }
  return ResOk;
}

void SamplesGeneratorImpl::int2Comb(std::vector<int> &oMaterials, int idx, int nMat, int nVox)
{
  oMaterials.resize(nVox);

  for(int ii = nVox-1;ii>=0;ii--){
    oMaterials[ii] = idx % nMat;
    idx/=nMat;
  }
}

int getVoxelIndex(int i, int j, int k, int nx, int ny, int nz)
{
  return i*ny*nz + j*nz + k;
}

void SamplesGeneratorImpl::int2Comb_biMaterialStructures(std::vector<int> &oMaterials, int idx, int nMat, int iN[3])
{
  int nMat2 = nMat*nMat;
  if (idx>=nMat2)
    idx++;
  if (idx>=2*nMat2)
    idx++;

  int indMat1 = (idx%nMat2)/nMat;
  int indMat2 = (idx%nMat2)%nMat;

   int n[3] = {2,2,2}; // {iN[0], iN[1], iN[2]};

  oMaterials.clear();
  oMaterials.resize(n[0]*n[1]*n[2], indMat2);

  int axis = idx / (nMat*nMat);
  n[axis] /= 2;

  int i,j,k;
  for (i=0; i<n[0]; i++)
  {
    for (j=0; j<n[1]; j++)
    {
      for (k=0; k<n[2]; k++)
      {
        int voxelIndex = getVoxelIndex(i,j,k, n[0], n[1], n[2]);
        oMaterials[voxelIndex] = indMat1;
      }
    }
  }
}

void SamplesGeneratorImpl::int2CombRandom(std::vector<int> &oMaterials,  int nMat, int nVox)
{
  oMaterials.resize(nVox);
  for(int ii=0; ii<nVox; ii++)
  {
    oMaterials[ii] = rand() % nMat;
  }
}

void SamplesGeneratorImpl::vectorMat2MatrixMat(const std::vector<int> &iVecMaterials, int oMatrixMaterials[2][2][2])
{
  for (int i=0; i<8; i++)
  {
    int indx = i/4;
    int indy = (i%4)/2;
    int indz = i%2;
    oMatrixMaterials[indx][indy][indz] = iVecMaterials[i];
  }
}

int SamplesGeneratorImpl::computeMaterialParameters(std::string iStepperType , const std::vector<std::vector<int> > &iMaterialAssignments, 
                                                    const std::vector<std::vector<int> > &iBaseMaterialStructures, int iBaseMatStructureSize[3], int iLevel, int iBlockRep)
{
  int c = 2*iBlockRep;
  int n[3] = {c*iBaseMatStructureSize[0], c*iBaseMatStructureSize[1], c*iBaseMatStructureSize[2]};
  int N[3] = {iBaseMatStructureSize[0], iBaseMatStructureSize[1], iBaseMatStructureSize[2]};

  int icomb, ncomb=(int)iMaterialAssignments.size();
  for (icomb=0; icomb<ncomb; icomb++)
  {
    std::vector<int> cellMaterials;
    if (m_dim==3)
    {
      getMaterialAssignment(n[0], n[1], n[2], iMaterialAssignments[icomb], N[0], N[1], N[2], iBaseMaterialStructures, cellMaterials);
    }
    else
    {
       getMaterialAssignment(n[0], n[1], iMaterialAssignments[icomb], N[0], N[1], iBaseMaterialStructures, cellMaterials);
    }
    int iaxis=0;
    for (iaxis=0; iaxis<2; iaxis++)
    {
      //float forceMagnitude = 1.f;
      float forceMagnitude = 0.3f;
      int NumberOfSample = m_UseLinearMaterial? 1: 10;
      int axis = iaxis;
      //int axis = 0;
      Vector3f forceAxis;
      forceAxis[axis] = 1.f;

      std::string strAxis =  (iaxis==0?"x//": "y//");
      std::string FileName = m_OutputDirectory + strAxis  + "StressStrain_" + std::to_string(iLevel) + '_' + std::to_string(3) + "_" + std::to_string(icomb) + ".txt";

      std::vector<float> stresses, strains;
      bool ResOk = true;
      if (m_dim==3)
      {
        sampleDeformation(n, m_mat, iStepperType, forceMagnitude, axis, NumberOfSample, cellMaterials, stresses, strains);
      }
      else
      {
        sampleDeformation(n, m_mat2D, iStepperType, forceMagnitude, axis, NumberOfSample, cellMaterials, stresses, strains);
      }
      if (ResOk)
      {
        writeData(FileName, forceAxis, iMaterialAssignments[icomb], stresses, strains);
      }
    }
  }
  //exit(0);
  return 0;
}

int SamplesGeneratorImpl::computeMaterialParametersLevel1(std::string & iStepperType)
{
  // level 0
  std::vector<std::vector<int> > baseMaterialStructures;
  baseMaterialStructures.resize(2);
  baseMaterialStructures[0].push_back(0);
  baseMaterialStructures[1].push_back(1);
  int N[3] = {1,1,1};

  // level 1
  int level = 1;
  int nVoxCoarse = pow(2, m_dim);
  int nmat = 2;
  int ncomb=pow(nmat,nVoxCoarse);
  std::vector<std::vector<int> > materialAssignments(ncomb);
  int icomb;
  for (icomb=0; icomb<ncomb; icomb++)
  {
    int2Comb(materialAssignments[icomb], icomb, nmat, nVoxCoarse);
  }

  std::string FileNameMaterialsBaseStructure = m_OutputDirectory + "Materials_" + std::to_string(level)  + ".txt";
  bool ResOk = writeMaterialCombinations(FileNameMaterialsBaseStructure,  baseMaterialStructures);

  computeMaterialParameters(iStepperType, materialAssignments, baseMaterialStructures, N, level, m_blockRep);

  return 0;
}

int SamplesGeneratorImpl::computeMaterialParameters(std::string & iStepperType, int iLevel, int iNbCombPrevLevel)
{
  assert(iLevel>=2);

  int prevLevel = iLevel-1;
  int blockSize = m_blockRep;

  std::string subDirectories[2];
  subDirectories[0] =  m_OutputDirectory + "x_level" + std::to_string(prevLevel) + "//";
  subDirectories[1] =  m_OutputDirectory + "y_level" + std::to_string(prevLevel) + "//";
  
  std::string baseFileName = "StressStrain_" + std::to_string(prevLevel) + "_" + std::to_string(blockSize);
  int NbFiles = iNbCombPrevLevel;

  std::vector<cfgScalar> physicalParameters; //YoungModulus, density;
  std::vector<std::vector<int> > materialAssignments, materialAssignmentsOneCell, baseMaterialStructures;

  std::string materialFile = m_OutputDirectory + "Materials_" + std::to_string(prevLevel) + ".txt";
  bool ResOk = cfgMaterialUtilities::computeMaterialParameters(m_dim, materialFile, subDirectories, baseFileName, NbFiles, physicalParameters, baseMaterialStructures, materialAssignments, materialAssignmentsOneCell);

  int naxis = 2;
  int nparam = naxis + 1;
  std::vector<int> convexHull, boundaryVertices, convexHullFaces;
  std::vector<cfgScalar> parameterPoints;
  parameterPoints.insert(parameterPoints.end(), physicalParameters.begin(), physicalParameters.end());
  //parameterPoints.insert(parameterPoints.end(), physicalParameters.begin(), physicalParameters.begin()+ 3*m_lastPointIndices[1]+3);
  computeDelaundayTriangulation(parameterPoints, nparam, convexHull, boundaryVertices, convexHullFaces);

  //computeConvexHull(physicalParameters, nparam, convexHullVertices);
  if (1)
  {
    writeVector2File(boundaryVertices, m_OutputDirectory + "convexHull.m");
  }
  std::vector<std::vector<int> > newBaseMaterialStructures;

  int sideSize = pow(2, prevLevel);
  int N[3] = {sideSize,sideSize,sideSize};
  std::vector<std::vector<int> > newMaterialAssignments;

  int version=2;
  if (version==0)
  {
    int nVoxCoarse = 8;
    int nmat = (int)boundaryVertices.size();
    //int nNewComb=pow(nmat,nVoxCoarse);
    int nNewComb = 3*nmat*nmat-2;
    newMaterialAssignments.resize(nNewComb);
    int icomb=0;
    for (icomb=0; icomb<nNewComb; icomb++)
    {
      //int2Comb(newMaterialAssignments[icomb], icomb, nmat, nVoxCoarse);
      int2Comb_biMaterialStructures(newMaterialAssignments[icomb], icomb, nmat, N);
    }
    newBaseMaterialStructures = getSubVector(materialAssignments, boundaryVertices);
    int borderSimilarity = 0;
    float avgBorderSimilarity = 0;
    for (int icomb=0; icomb<nNewComb; icomb++)
    {
      borderSimilarity = checkBorderSimilarity(newMaterialAssignments[icomb], N,  newBaseMaterialStructures);
      avgBorderSimilarity += borderSimilarity;
    }
    avgBorderSimilarity /= (float)nNewComb;
  }
  else if (version==1)
  {
    int nNewComb = 1000;
    std::vector<int> BaseMaterialStructuresIndices;
    getNewCombinations(boundaryVertices, physicalParameters, N, materialAssignments, nNewComb, newMaterialAssignments, BaseMaterialStructuresIndices);
    //newBaseMaterialStructures = getSubVector(materialAssignments, BaseMaterialStructuresIndices);
    newBaseMaterialStructures = materialAssignments;

    int borderSimilarity = 0;
    float avgBorderSimilarity = 0;
    for (int icomb=0; icomb<nNewComb; icomb++)
    {
      borderSimilarity = checkBorderSimilarity(newMaterialAssignments[icomb], N,  newBaseMaterialStructures);
      avgBorderSimilarity += borderSimilarity;
    }
    avgBorderSimilarity /= (float)nNewComb;
  }
  else if (version==2) //exhaustive search
  {
    int nmat = (int)materialAssignments.size();
    int nVoxCoarse = pow(2, m_dim);
    int nNewComb=pow(nmat,nVoxCoarse);
    newMaterialAssignments.resize(nNewComb);
    int icomb;
    for (icomb=0; icomb<nNewComb; icomb++)
    {
      int2Comb(newMaterialAssignments[icomb], icomb, nmat, nVoxCoarse);
    }
    newBaseMaterialStructures = materialAssignments;
  }
  else if (version==3)
  {
    /*int nNewComb = 1000;
    std::vector<int> BaseMaterialStructuresIndices;
    getNewCombinationsV2(boundaryVertices, physicalParameters, N, materialAssignments, nNewComb, newMaterialAssignments, BaseMaterialStructuresIndices);
    newBaseMaterialStructures = materialAssignments;*/ 
  }

  std::string FileNameMaterials = m_OutputDirectory + "Materials_" + std::to_string(iLevel)  + ".txt";
  writeMaterialCombinations(FileNameMaterials,  newBaseMaterialStructures);

  int blockRep = 3;
  computeMaterialParameters(iStepperType, newMaterialAssignments, newBaseMaterialStructures, N, iLevel, blockRep);

  return 0;
}

void SamplesGeneratorImpl::getNewCombinations(const std::vector<int> &iBoundaryPoints, const std::vector<cfgScalar> &iPoints, 
                                              int N[3], const std::vector<std::vector<int> > &iBaseMaterialStructures, int iNbCombinations, 
                                              std::vector<std::vector<int> > &oNewMaterialAssignments, std::vector<int> &oNewBaseMaterialStructures)
{
  std::vector<cfgScalar> boundaryPointCoords = getSubVector(iPoints, 3, iBoundaryPoints);
  cfgScalar bary[3];
  getBarycenter<cfgScalar,3>(boundaryPointCoords, bary);
  Vector3f c(bary[0],bary[1],bary[2]);

  cfgScalar SqMaxDist = 0;
  int ipoint=0, npoint=(int)iBoundaryPoints.size();
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    int indPoint = iBoundaryPoints[ipoint];
    Vector3f p = getVector3f(indPoint, iPoints);
    cfgScalar SqDist = (p-c).absSquared();
    if (SqDist>SqMaxDist)
    {
      SqMaxDist = SqDist;
    }
  }

  cfgScalar SqMaxDist2 = 0;
  npoint = iPoints.size()/3;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    Vector3f p = getVector3f(ipoint, iPoints);
    int closestBoundaryPoint = getClosestPoint(p, boundaryPointCoords);
    Vector3f q = getVector3f(iBoundaryPoints[closestBoundaryPoint], iPoints);
    float SqDist = (p-q).absSquared();
    if (SqDist>SqMaxDist2)
    {
      SqMaxDist2 = SqDist;
    }
  }
  SqMaxDist = SqMaxDist2;

  int maxBoundaryMismatch = 3*N[0]*N[0]*4*2;
 
  oNewBaseMaterialStructures.clear();
  oNewMaterialAssignments.clear();
  int n[3] = {2,2,2};
  int nMat = iPoints.size()/3;
  int nVox = n[0]*n[1]*n[2];
  std::vector<int> materialAssignement(nVox, -1);
  float pmin = 0.1;
  float pmax = 1;
  std::set<int> BaseMaterials;
  std::vector<int> voxelIndices;
  int ivox;
  for (ivox=0; ivox<nVox; ivox++)
  {
    voxelIndices.push_back(ivox);
  }
  while (oNewMaterialAssignments.size()<iNbCombinations)
  {
    materialAssignement.clear();
    materialAssignement.resize(nVox, -1);

    float r2 = (float)rand()/float(RAND_MAX);
    random_shuffle(voxelIndices.begin(), voxelIndices.end());

    bool assignmentOk = true;
    int ivox=0;
    for (ivox=0; ivox<nVox && assignmentOk; ivox++)
    {
      int indVox = voxelIndices[ivox];

      bool matAssigned = false;  
      while (!matAssigned && assignmentOk)
      {
        int matIndex = rand() % nMat;
        Vector3f p = getVector3f(matIndex, iPoints);
        int closestBoundaryPoint = getClosestPoint(p, boundaryPointCoords);
        Vector3f q = getVector3f(iBoundaryPoints[closestBoundaryPoint], iPoints);
        float SqDist = (p-q).absSquared();
        float a = sqrt(SqDist/SqMaxDist);
        float proba = a*pmin + (1-a)*pmax;
        /*if (SqDist < 1.e-6)
        {
        proba = 1;
        }
        else
        {
        proba = 0.1;
        }*/ 
        materialAssignement[indVox] = matIndex;
        int borderSimilarity = checkBorderSimilarity(materialAssignement, N,  iBaseMaterialStructures);
        float ratioMismatch = (float)borderSimilarity/(float)maxBoundaryMismatch;
        float proba2 = 1-ratioMismatch;
        //proba2 *= proba2*proba2;
        if (ratioMismatch > 0.34)
        {
          proba2 = 0.1;
        }
        else
        {
          proba2 = proba2;
        }

        float r = (float)rand()/float(RAND_MAX);

        if (proba2 < r2)
        {
          assignmentOk = false;
        }
        else if (proba > r)
        {
          materialAssignement[indVox] = matIndex;
          matAssigned = true;
          BaseMaterials.insert(matIndex);
        }
      }
    }
    if (assignmentOk)
    {
      oNewMaterialAssignments.push_back(materialAssignement);
      std::cout << "# assignments = " << oNewMaterialAssignments.size() << std::endl;
    }
  }
  oNewBaseMaterialStructures = toStdVector(BaseMaterials);
}

int SamplesGeneratorImpl::getClosestPoint(Vector3f &iP, const std::vector<cfgScalar> &iPoints)
{
  int ClosestPointIndex = -1;
  float SqDistMin = FLT_MAX;
  int ipoint=0, npoint=(int)iPoints.size()/3;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    Vector3f Q = getVector3f(ipoint, iPoints);
    float SqDist = (Q-iP).absSquared();
    if (SqDist < SqDistMin)
    {
      SqDistMin = SqDist;
      ClosestPointIndex = ipoint;
    }
  }
  return ClosestPointIndex;
}

int SamplesGeneratorImpl::run()
{
  srand(0);

  //std::vector<StrainEneNeo> ene(2);
  std::vector<StrainLin2D> ene2D(2);
  ene2D[0].param[0] = 1;
  ene2D[0].param[1] = 10;
  ene2D[1].param[0] = 10;
  ene2D[1].param[1] = 100;

  //std::vector<StrainEneNeo> ene(2);
  std::vector<StrainLin> ene(2);
  ene[0].param[0] = 1;
  ene[0].param[1] = 10;
  ene[1].param[0] = 10;
  ene[1].param[1] = 100;

  if (m_dim==2)
  {
    m_mat2D.clear();
    m_mat2D.resize(ene2D.size());
    for (unsigned int ii = 0; ii < m_mat2D.size(); ii++){
      for (unsigned int jj = 0; jj < m_mat2D[ii].e.size(); jj++){
        m_mat2D[ii].e[jj] = &ene2D[ii];
      }
    }
  }
  else
  {
   
    m_mat.clear();
    m_mat.resize(ene.size());
    for (unsigned int ii = 0; ii < m_mat.size(); ii++){
      for (unsigned int jj = 0; jj < m_mat[ii].e.size(); jj++){
        m_mat[ii].e[jj] = &ene[ii];
      }
    }
  }
  //std::string stepperType = "newtonCuda";
  std::string stepperType = "newton";

  if (1)
  {
    //computeMaterialParametersLevel1(stepperType);
    //computeMaterialParameters(stepperType, 2, (m_dim==2?16:256));
    computeMaterialParameters(stepperType, 3, 20);

    //computeMaterialParameters(material, stepperType, 3, 586);
    //computeMaterialParameters(material, stepperType, 3, 429);

    return 0 ;
  }

  //sampleMaterialSpace();


  /*std::vector<float> points;
  points.push_back(-1);
  points.push_back(-1);

  points.push_back(-1);
  points.push_back(1);

  points.push_back(1);
  points.push_back(-1);

  points.push_back(1);
  points.push_back(1);
 
  computeConvexHull(points, 2);*/ 

  int shiftRand = 0;
  for (int i=0; i<shiftRand; i++)
  {
    rand();
  }
 
  int MinHierarchyLevel = 2;
  int MaxHierarchyLevel = 2;

  std::vector<std::vector<int> > materialCombinations;
  if (MinHierarchyLevel<=1)
  {
    materialCombinations.resize(2);
    materialCombinations[0].push_back(0);
    materialCombinations[1].push_back(1);
  }
  else if (MinHierarchyLevel==2)
  {
    int nVoxCoarse = 8;
    int nmat = 2;
    int icomb=0, ncomb=pow(nmat,nVoxCoarse);
     materialCombinations.resize(ncomb);
    for (icomb=0; icomb<ncomb; icomb++)
    {
      int2Comb(materialCombinations[icomb], icomb, nmat, nVoxCoarse);
    }
  }

  /*int iHLevel;
  for (iHLevel=MinHierarchyLevel; iHLevel<=MaxHierarchyLevel; iHLevel++)
  {
    std::string FileNameMaterials =  m_OutputDirectory + "Materials_" + std::to_string(iHLevel)  + ".txt";
    writeMaterialCombinations(FileNameMaterials,  materialCombinations);

    int nmat = materialCombinations.size();
    int nVoxCoarse = 8;
    int nCoarse = pow(2, iHLevel-1);

    int nMaxComb = 400;
    int icomb=0;
    int ncomb = (iHLevel<2? ncomb=pow(nmat,nVoxCoarse): nMaxComb);
    ncomb=std::min(nMaxComb, ncomb);
    std::vector<std::vector<int> > newmaterialCombinations(ncomb);
    for (icomb=0; icomb<ncomb; icomb++)
    {
      std::vector<int> coarseCellMaterials;
      int2Comb(coarseCellMaterials, icomb, nmat, nVoxCoarse);
      //int2CombRandom(coarseCellMaterials, icomb, nmat, nVoxCoarse);
      //coarseCellMaterials[0] = 1;

      //vectorMat2MatrixMat(mat, macroStructureMaterials);

      int MinRefinementLevel = iHLevel+1;
      int MaxRefinementLevel = iHLevel+1;
      int iReflevel;
      for (iReflevel=MinRefinementLevel; iReflevel<=MaxRefinementLevel; iReflevel++)
      {
        int res = (int)std::pow(2, iReflevel);
        int n[3] = {res, res, res};
        int nVox = n[0]*n[1]*n[2];

        std::vector<int> gridMaterials;
        getMaterialAssignment(n[0], n[1], n[2], coarseCellMaterials, nCoarse, nCoarse, nCoarse, materialCombinations, gridMaterials);
        if (iReflevel==MinRefinementLevel)
        {
          newmaterialCombinations[icomb] = gridMaterials;
        }

        //float forceMagnitude = 1.f;
        float forceMagnitude = 0.3f;
        int NumberOfSample = 10;
        int axis = 0;
        Vector3f forceAxis;
        forceAxis[axis] = 1.f;

        std::string FileName = m_OutputDirectory + "StressStrain_" + std::to_string(iHLevel) + '_' + std::to_string(iReflevel) + "_" + std::to_string(icomb+shiftRand) + ".txt";

        std::vector<float> stresses, strains;
        sampleDeformation(n, material, stepperType, forceMagnitude, axis, NumberOfSample, gridMaterials, stresses, strains);
        writeData(FileName, forceAxis, coarseCellMaterials, stresses, strains);
      }
    }
    materialCombinations = newmaterialCombinations;
  }*/ 

  /*World * world = new World();
  world->em.push_back(em);

  Render render;
  render.init(world);
  render.loop();*/ 
  return 0;
}




