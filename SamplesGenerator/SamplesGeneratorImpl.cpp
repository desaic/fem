
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
#include "StepperGrad.hpp"

#include "MeshUtilities.h"
using namespace meshUtil;

#include <algorithm>

#include "cfgUtilities.h"
using namespace cfgUtil;

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

#include "ScoringFunction.h"
#include "Resampler.h"
#include "NumericalCoarsening.h"
#include "MaterialOptimizer.h"
#include "DistanceField.h"

SamplesGeneratorImpl::SamplesGeneratorImpl(int iDim)
{
  m_UseLinearMaterial = true;
  assert(iDim==2||iDim==3);
  //m_dim = iDim;
  m_dim = 3;
  m_blockRep = 1;// (m_dim==2? 8: 4);
  m_nbSubdivisions = 1;
  m_orthotropicOnly = false;
  m_cubicOnly = true;
  m_continuousMatDist = true;
}

SamplesGeneratorImpl::~SamplesGeneratorImpl()
{
}

float SamplesGeneratorImpl::computeStrain(const ElementRegGrid * iElementGrid, int iAxis)
{
  //writeVector2File(iElementGrid->x, m_OutputDirectory + "x.m");

  Vector3S Barycenters[2];

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
  float strain = (Barycenters[0]-Barycenters[1]).norm()-1;
  return strain;
}

float SamplesGeneratorImpl::computeStrain(const ElementRegGrid2D * iElementGrid, int iAxis)
{
  //writeVector2File(iElementGrid->x, m_OutputDirectory + "x.m");

  Vector2S Barycenters[2];

  int iside;
  for (iside=0; iside<2; iside++)
  {
    Barycenters[iside ] = Vector2S::Zero();
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
  float strain = (Barycenters[0]-Barycenters[1]).norm()-1;
  return strain;
}


void SamplesGeneratorImpl::setExternalForces(ElementRegGrid * iElementGrid, int iAxis, int iSide, float iForceMagnitude)
{
  assert(iElementGrid);

  ElementRegGrid * em = iElementGrid;
 
  int n[3];
  n[0] = em->nx;
  n[1] = em->ny;
  n[2] = em->nz;

  Vector3S ff(0,0,0);
  ff[iAxis] = (1.f/(n[(iAxis+1)%3]*n[(iAxis+2)%3])) * iForceMagnitude/4;

  std::vector<int> sides;
  std::vector<int> signs;
  if (iSide==0)
  {
    sides.push_back(0);
    signs.push_back(-1);
  }
  else if (iSide==1)
  {
    sides.push_back(1);
    signs.push_back(1);
  }
  else
  {
    sides.push_back(0);
    sides.push_back(1);

    signs.push_back(-1);
    signs.push_back(1);
  }

  int iside=0, nside=(int)signs.size();
  for (iside=0; iside<nside; iside++)
  {
    int indSide = sides[iside];
    int sign = signs[iside];
    std::vector<int> elemIndices;
    std::vector<std::vector<int> > fvIndices;
    getSideVertices(2*iAxis+indSide, em, elemIndices, fvIndices);

    int ielem=0, nelem=(int)elemIndices.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      int indElement = elemIndices[ielem];
      int ivertex=0, nvertex=(int)fvIndices[ielem].size();
      for (ivertex=0; ivertex<nvertex; ivertex++)
      {
        int FvIndex = fvIndices[ielem][ivertex];
        int VertexIndex =em->e[indElement]->at(FvIndex);
        em->fe[VertexIndex] += sign*ff;
      }
    }
  }
  
  // Constraints
  if (nside==1)
  {
    std::vector<int> elemIndicesOppositeSide;
    std::vector<std::vector<int> > fvIndicesOppositeSide;
    getSideVertices(2*iAxis+1-iSide, em, elemIndicesOppositeSide, fvIndicesOppositeSide);
    int nelemOppositeSide=(int)elemIndicesOppositeSide.size();
    int ielem = 0;
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
}

void SamplesGeneratorImpl::setExternalForces(ElementRegGrid2D * iElementGrid, int iAxis, int iSide, float iForceMagnitude)
{
  assert(iElementGrid);

  ElementRegGrid2D * em = iElementGrid;

  int n[2];
  n[0] = em->nx;
  n[1] = em->ny;

  double thickness = 1.;

  Vector2S ff(0,0);
  ff[iAxis] = (1.f/(n[(iAxis+1)%2]*thickness)) * iForceMagnitude/2;

  std::vector<int> sides;
  std::vector<int> signs;
  if (iSide==0)
  {
    sides.push_back(0);
    signs.push_back(-1);
  }
  else if (iSide==1)
  {
    sides.push_back(1);
    signs.push_back(1);
  }
  else
  {
    sides.push_back(0);
    sides.push_back(1);

    signs.push_back(-1);
    signs.push_back(1);
  }

  int iside=0, nside=(int)signs.size();
  for (iside=0; iside<nside; iside++)
  {
    int indSide = sides[iside];
    int sign = signs[iside];
    std::vector<int> elemIndices;
    std::vector<std::vector<int> > fvIndices;
    getSideVertices(2*iAxis+indSide, em, elemIndices, fvIndices);

    int ielem=0, nelem=(int)elemIndices.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      int indElement = elemIndices[ielem];
      int ivertex=0, nvertex=(int)fvIndices[ielem].size();
      for (ivertex=0; ivertex<nvertex; ivertex++)
      {
        int FvIndex = fvIndices[ielem][ivertex];
        int VertexIndex =em->e[indElement]->at(FvIndex);
        em->fe[VertexIndex] += sign*ff;
      }
    }
  }

  // Constraints
  if (nside==1)
  {
    std::vector<int> elemIndicesOppositeSide;
    std::vector<std::vector<int> > fvIndicesOppositeSide;
    getSideVertices(2*iAxis+1-iSide, em, elemIndicesOppositeSide, fvIndicesOppositeSide);
    int nelemOppositeSide=(int)elemIndicesOppositeSide.size();
    int ielem = 0;
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
    ((StepperNewton2D*)stepper)->removeTranslation(true);
    ((StepperNewton2D*)stepper)->removeRotation(true);
    ((StepperNewton2D*)stepper)->enforcePeriodicity(true);
  }
  else if (iStepperType == "newtonCuda"){
 //   stepper = new NewtonCuda();
 //   ((NewtonCuda*)stepper)->m_fixRigidMotion = false;
  }
  else if (iStepperType == "ipopt"){
    //not yet implemented
    //  stepper = new IpoptStepper();
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
    ((StepperNewton2D*)stepper)->removeTranslation(true);
    ((StepperNewton2D*)stepper)->removeRotation(true);
    ((StepperNewton2D*)stepper)->enforcePeriodicity(true);
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
     getSideVoxels(2*iaxis, N[0], N[1], N[2], voxels1);
     getSideVoxels(2*iaxis+1, N[0], N[1], N[2], voxels2);

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
    Vector3S forceAxis;
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

  std::vector<Vector3S> xinit;
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
      //std::cout << "Strain = " << strain << std::endl;

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

  std::vector<Vector2S> xinit;
  int isample;
  for (isample=0; isample<iNumberOfSample && ResOk; isample++)
  {
    float forceMagnitude = (isample+1.f)*iMaxForce/iNumberOfSample;

    ElementRegGrid2D * em = createPhysicalSystem(iN, iMaterial);
    setExternalForces(em, iForceAxis, 2, forceMagnitude);
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
      //std::cout << "Strain = " << strain << std::endl;

      delete em;

      oStresses.push_back(forceMagnitude);
      oStrains.push_back(strain);
    }
  }
  return ResOk;
}

bool SamplesGeneratorImpl::sampleDeformation(int iN[2], std::vector<MaterialQuad2D> &iMaterial, const std::string iStepperType, float iMaxForce, int iForceAxis, int iNumberOfSample, const std::vector<int> & iMaterials,
                                             std::vector<float> &oStresses, std::vector<std::vector<float> > &oDeformations)
{
  bool ResOk = true;

  oStresses.clear();
  oDeformations.clear();

  std::vector<Vector2S> xinit;
  int isample;
  for (isample=0; isample<iNumberOfSample && ResOk; isample++)
  {
    float forceMagnitude = (isample+1.f)*iMaxForce/iNumberOfSample;

    ElementRegGrid2D * em = createPhysicalSystem(iN, iMaterial);
    setExternalForces(em, iForceAxis, 2, forceMagnitude);
    em->me = iMaterials;
    //assignMaterials(em, iMaterials);

    if (isample>0)
    {
      em->x = xinit;
    }
    bool Converged = false;
    ResOk = computeForceDeformationSample(em, iStepperType, 1, Converged);
    if (ResOk)
    {
      xinit = em->x;

      if (0)
      {
        float strain = computeStrain(em, iForceAxis);
      }

      //writeVector2File(em->me, m_OutputDirectory + "materials.m");
      //std::cout << "Strain = " << strain << std::endl;

      delete em;

      oStresses.push_back(forceMagnitude);
      oDeformations.push_back(toVectorScalar(xinit));
    }
  }
  return ResOk;
}

bool SamplesGeneratorImpl::sampleDeformation(int iN[3], std::vector<MaterialQuad> &iMaterial, const std::string iStepperType, float iMaxForce, int iForceAxis, int iNumberOfSample, const std::vector<int> & iMaterials,
                                             std::vector<float> &oStresses, std::vector<std::vector<float> > &oDeformations)
{
  bool ResOk = true;

  oStresses.clear();
  oDeformations.clear();

  std::vector<Vector3S> xinit;
  int isample;
  for (isample=0; isample<iNumberOfSample && ResOk; isample++)
  {
    float forceMagnitude = (isample+1.f)*iMaxForce/iNumberOfSample;

    ElementRegGrid * em = createPhysicalSystem(iN, iMaterial);
    setExternalForces(em, iForceAxis, 2, forceMagnitude);
    em->me = iMaterials;
    //assignMaterials(em, iMaterials);

    if (isample>0)
    {
      em->x = xinit;
    }
    bool Converged = false;
    ResOk = computeForceDeformationSample(em, iStepperType, 1, Converged);
    if (ResOk)
    {
      xinit = em->x;

      if (0)
      {
        float strain = computeStrain(em, iForceAxis);
      }

      //writeVector2File(em->me, m_OutputDirectory + "materials.m");
      //std::cout << "Strain = " << strain << std::endl;

      delete em;

      oStresses.push_back(forceMagnitude);
      oDeformations.push_back(toVectorFloat(xinit));
    }
  }
  return ResOk;
}

void SamplesGeneratorImpl::int2Comb(std::vector<int> &oMaterials, int idx, int nMat, int nVox)
{
  oMaterials.resize(nVox);

  //for(int ii = nVox-1;ii>=0;ii--){
  for(int ii = 0; ii<nVox; ii++)
  {
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
                                                    const std::vector<std::vector<int> > &iBaseMaterialStructures, int iBaseMatStructureSize[3], int iLevel, int iBlockRep, bool iWriteSingleFile, bool iWriteFullDeformation)
{
  int c = 2*iBlockRep*m_nbSubdivisions;
  int n[3] = {c*iBaseMatStructureSize[0], c*iBaseMatStructureSize[1], m_dim==2? 0: c*iBaseMatStructureSize[2]};
  int N[3] = {iBaseMatStructureSize[0], iBaseMatStructureSize[1], iBaseMatStructureSize[2]};

  int ncomb=(int)iMaterialAssignments.size();
  
  std::vector<std::vector<float> > stresses[2], strains[2];
  std::vector<std::vector<std::vector<float> > > x[2];
  int iaxis=0;
  for (iaxis=0; iaxis<2; iaxis++)
  {
    stresses[iaxis].resize(ncomb);
    strains[iaxis].resize(ncomb);
    x[iaxis].resize(ncomb);
  }

  for (int icomb=0; icomb<ncomb; icomb++)
  {
    std::vector<int> cellMaterials;
    if (m_dim==3)
    {
      getMaterialAssignment(2, 2, 2, iMaterialAssignments[icomb], N[0], N[1], N[2], iBaseMaterialStructures, iBlockRep, iBlockRep, iBlockRep, m_nbSubdivisions, cellMaterials);
    }
    else
    {
       getMaterialAssignment(2, 2, iMaterialAssignments[icomb], N[0], N[1], iBaseMaterialStructures, iBlockRep, iBlockRep, m_nbSubdivisions, cellMaterials);
    }
    int iaxis=0;
    for (iaxis=0; iaxis<2; iaxis++)
    {
      float forceMagnitude = 1.f;
      int NumberOfSample = m_UseLinearMaterial? 1: 10;
      int axis = iaxis;
  
      if (m_dim==3)
      {
        sampleDeformation(n, m_mat, iStepperType, forceMagnitude, axis, NumberOfSample, cellMaterials, stresses[iaxis][icomb], strains[iaxis][icomb]);
      }
      else
      {
        if (iWriteFullDeformation)
        {
          sampleDeformation(n, m_mat2D, iStepperType, forceMagnitude, axis, NumberOfSample, cellMaterials, stresses[iaxis][icomb], x[iaxis][icomb]);
        }
        else
        {
          sampleDeformation(n, m_mat2D, iStepperType, forceMagnitude, axis, NumberOfSample, cellMaterials, stresses[iaxis][icomb], strains[iaxis][icomb]);
        }
      }

      if (icomb%1000==0)
      {
        std::cout << "ncomb = " << icomb << std::endl;
        Vector3S forceAxis(0,0,0);
        forceAxis[iaxis] = 1.f;

        std::vector<std::vector<int> > materials = iMaterialAssignments;
        materials.resize(icomb);

        if (iWriteFullDeformation)
        {
          std::string FileName = m_OutputDirectory + "StressDeformation_" + std::to_string(iLevel) + '_' + std::to_string(m_blockRep) + "_" + (iaxis==0? "x": "y") + ".txt";
          //writeData(FileName, forceAxis, materials, stresses[iaxis], x[iaxis], n);
          writeDataBinary(FileName, forceAxis, materials, stresses[iaxis], x[iaxis], n);

        }
        else
        {
          std::string FileName = m_OutputDirectory + "StressStrain_" + std::to_string(iLevel) + '_' + std::to_string(m_blockRep) + "_" + (iaxis==0? "x": "y") + ".txt";
          writeData(FileName, forceAxis, materials, stresses[iaxis], strains[iaxis]);
        }
      }
    }
  }
  if (iWriteSingleFile)
  {
    int iaxis=0;
    for (iaxis=0; iaxis<2; iaxis++)
    {
      Vector3S forceAxis(0,0,0);
      forceAxis[iaxis] = 1.f;
      if (iWriteFullDeformation)
      {
        std::string FileName = m_OutputDirectory + "StressDeformation_" + std::to_string(iLevel) + '_' + std::to_string(m_blockRep) + "_" + (iaxis==0? "x": "y") + ".txt";
        writeData(FileName, forceAxis, iMaterialAssignments, stresses[iaxis], x[iaxis], n);
      }
      else
      {
        std::string FileName = m_OutputDirectory + "StressStrain_" + std::to_string(iLevel) + '_' + std::to_string(m_blockRep) + "_" + (iaxis==0? "x": "y") + ".txt";
        writeData(FileName, forceAxis, iMaterialAssignments, stresses[iaxis], strains[iaxis]);
      }
    }
  }
  else
  {
    for (int icomb=0; icomb<ncomb; icomb++)
    {
      int iaxis=0;
      for (iaxis=0; iaxis<2; iaxis++)
      {
        Vector3S forceAxis(0,0,0);
        forceAxis[iaxis] = 1.f;
        std::string strAxis =  (iaxis==0?"x//": "y//");
        std::string FileName = m_OutputDirectory + strAxis  + "StressStrain_" + std::to_string(iLevel) + '_' + std::to_string(3) + "_" + std::to_string(icomb) + ".txt";
        writeData(FileName, forceAxis, iMaterialAssignments[icomb], stresses[iaxis][icomb], strains[iaxis][icomb]);
      }
    }
  }
  //exit(0);
  return 0;
}

int SamplesGeneratorImpl::computeMaterialParameters(std::string iStepperType , const std::vector<std::vector<int> > &iMaterialAssignments, int iNewMatStructureSize[3],
                                                    const std::vector<std::vector<int> > &iBaseMaterialStructures, int iBaseMatStructureSize[3], int iLevel, int iBlockRep, bool iWriteSingleFile, bool iWriteFullDeformation,
                                                    const std::string &iPostfix)
{
  int c = iBlockRep;
  int n[3] = {c*iNewMatStructureSize[0]*m_nbSubdivisions, c*iNewMatStructureSize[1]*m_nbSubdivisions, m_dim==2? 0: c*iNewMatStructureSize[2]*m_nbSubdivisions};
  int N[3] = {iBaseMatStructureSize[0], iBaseMatStructureSize[1], iBaseMatStructureSize[2]};

  int ncomb=(int)iMaterialAssignments.size();
  
  std::vector<std::vector<float> > stresses[2], strains[2];
  std::vector<std::vector<std::vector<float> > > x[2];
  int iaxis=0;
  for (iaxis=0; iaxis<2; iaxis++)
  {
    stresses[iaxis].resize(ncomb);
    strains[iaxis].resize(ncomb);
    x[iaxis].resize(ncomb);
  }

  for (int icomb=0; icomb<ncomb; icomb++)
  {
    std::vector<int> cellMaterials;
    if (m_dim==3)
    {
      getMaterialAssignment(iLevel, iLevel, iLevel, iMaterialAssignments[icomb], N[0], N[1], N[2], iBaseMaterialStructures, iBlockRep, iBlockRep, iBlockRep, m_nbSubdivisions, cellMaterials);
    }
    else
    {
       getMaterialAssignment(iLevel, iLevel, iMaterialAssignments[icomb], N[0], N[1], iBaseMaterialStructures, iBlockRep, iBlockRep, m_nbSubdivisions, cellMaterials);
    } 
    //cellMaterials = iMaterialAssignments[icomb];

    int iaxis=0;
    for (iaxis=0; iaxis<2; iaxis++)
    {
      float forceMagnitude = 1.f;
      int NumberOfSample = m_UseLinearMaterial? 1: 10;
      int axis = iaxis;
  
      if (m_dim==3)
      {
        sampleDeformation(n, m_mat, iStepperType, forceMagnitude, axis, NumberOfSample, cellMaterials, stresses[iaxis][icomb], strains[iaxis][icomb]);
      }
      else
      {
        if (iWriteFullDeformation)
        {
          sampleDeformation(n, m_mat2D, iStepperType, forceMagnitude, axis, NumberOfSample, cellMaterials, stresses[iaxis][icomb], x[iaxis][icomb]);
        }
        else
        {
          sampleDeformation(n, m_mat2D, iStepperType, forceMagnitude, axis, NumberOfSample, cellMaterials, stresses[iaxis][icomb], strains[iaxis][icomb]);
        }
      }
      if (icomb%50==0)
      {
        std::cout << "ncomb = " << icomb << std::endl;
        Vector3S forceAxis(0,0,0);
        forceAxis[iaxis] = 1.f;

        std::vector<std::vector<int> > materials = iMaterialAssignments;
        materials.resize(icomb);

        if (iWriteFullDeformation)
        {
          std::string extension = ".bin"; // ".txt";
          std::string FileName = m_OutputDirectory + "StressDeformation_" + std::to_string(iLevel) + '_' + std::to_string(m_blockRep) + "_" + (iaxis==0? "x": "y") + iPostfix + extension;
          //writeData(FileName, forceAxis, materials, stresses[iaxis], x[iaxis], n);
          writeDataBinary(FileName, forceAxis, materials, stresses[iaxis], x[iaxis], n);
        }
        else
        {
          std::string FileName = m_OutputDirectory + "StressStrain_" + std::to_string(iLevel) + '_' + std::to_string(m_blockRep) + "_" + (iaxis==0? "x": "y") + iPostfix + ".txt";
          writeData(FileName, forceAxis, materials, stresses[iaxis], strains[iaxis]);
        }
      }
    }
  }
  if (iWriteSingleFile)
  {
    int iaxis=0;
    for (iaxis=0; iaxis<2; iaxis++)
    {
      Vector3S forceAxis(0,0,0);
      forceAxis[iaxis] = 1.f;
      if (iWriteFullDeformation)
      {
        std::string extension = ".bin"; // ".txt";
        std::string FileName = m_OutputDirectory + "StressDeformation_" + std::to_string(iLevel) + '_' + std::to_string(m_blockRep) + "_" + (iaxis==0? "x": "y") + iPostfix + extension;
        //writeData(FileName, forceAxis, iMaterialAssignments, stresses[iaxis], x[iaxis], n);
        writeDataBinary(FileName, forceAxis, iMaterialAssignments, stresses[iaxis], x[iaxis], n);
      }
      else
      {
        std::string FileName = m_OutputDirectory + "StressStrain_" + std::to_string(iLevel) + '_' + std::to_string(m_blockRep) + "_" + (iaxis==0? "x": "y") + iPostfix + ".txt";
        writeData(FileName, forceAxis, iMaterialAssignments, stresses[iaxis], strains[iaxis]);
      }
    }
  }
  else
  {
    for (int icomb=0; icomb<ncomb; icomb++)
    {
      int iaxis=0;
      for (iaxis=0; iaxis<2; iaxis++)
      {
        Vector3S forceAxis(0,0,0);
        forceAxis[iaxis] = 1.f;
        std::string strAxis =  (iaxis==0?"x//": "y//");
        std::string FileName = m_OutputDirectory + strAxis  + "StressStrain_" + std::to_string(iLevel) + '_' + std::to_string(3) + "_" + std::to_string(icomb) + iPostfix + ".txt";
        writeData(FileName, forceAxis, iMaterialAssignments[icomb], stresses[iaxis][icomb], strains[iaxis][icomb]);
      }
    }
  }
  //exit(0);
  return 0;
}

int SamplesGeneratorImpl::computeDeformation(std::string iStepperType , const std::vector<std::vector<int> > &iMaterialAssignments, int iMatStructureSize[3],
                                             const std::vector<std::vector<int> > &iBaseMaterialStructures, int iBaseMatStructureSize[3], 
                                             std::vector<std::vector<float> > oStresses[3], std::vector<std::vector<std::vector<float> > > oX[3])
{
  int blockRep = m_blockRep;
  int n[3] = {blockRep*iMatStructureSize[0]*m_nbSubdivisions, blockRep*iMatStructureSize[1]*m_nbSubdivisions, m_dim==2? 0: blockRep*iMatStructureSize[2]*m_nbSubdivisions};
  int N[3] = {iBaseMatStructureSize[0], iBaseMatStructureSize[1], iBaseMatStructureSize[2]};

  int ncomb=(int)iMaterialAssignments.size();
  
  int iaxis=0;
  for (iaxis=0; iaxis<3; iaxis++)
  {
    oStresses[iaxis].resize(ncomb);
    oX[iaxis].resize(ncomb);
  }

  for (int icomb=0; icomb<ncomb; icomb++)
  {
    std::vector<int> cellMaterials;
    if (m_dim==3)
    {
      getMaterialAssignment(iMatStructureSize[0], iMatStructureSize[1], iMatStructureSize[2], iMaterialAssignments[icomb], N[0], N[1], N[2], iBaseMaterialStructures, blockRep, blockRep, blockRep, m_nbSubdivisions, cellMaterials);
    }
    else
    {
       getMaterialAssignment(iMatStructureSize[0], iMatStructureSize[1], iMaterialAssignments[icomb], N[0], N[1], iBaseMaterialStructures, blockRep, blockRep, m_nbSubdivisions, cellMaterials);
    } 
    int iaxis=0;
    for (iaxis=0; iaxis<m_dim; iaxis++)
    {
      float forceMagnitude = 1.f;
      int NumberOfSample = m_UseLinearMaterial? 1: 10;
      int axis = iaxis;
  
      if (m_dim==3)
      {
        sampleDeformation(n, m_mat, iStepperType, forceMagnitude, axis, NumberOfSample, cellMaterials, oStresses[iaxis][icomb], oX[iaxis][icomb]);
      }
      else
      {
        sampleDeformation(n, m_mat2D, iStepperType, forceMagnitude, axis, NumberOfSample, cellMaterials, oStresses[iaxis][icomb], oX[iaxis][icomb]);
      }
    }
    if (icomb%50==0)
    {
      std::cout << "ncomb = " << icomb << std::endl;
      std::vector<std::vector<int> > materials = iMaterialAssignments;
      materials.resize(icomb);

      std::vector<float> physicalParameters;
      if (m_dim==2)
      {
        std::vector<std::vector<float> > strains[2][2];
        cfgMaterialUtilities::computeStrain(n, oX, strains);
        cfgMaterialUtilities::computeMaterialParameters(iMaterialAssignments, iBaseMaterialStructures, oStresses, strains, physicalParameters);
      }
      else
      {
        std::vector<std::vector<float> > strains[3][3];
        cfgMaterialUtilities::computeStrain3D(n, oX, strains);
        cfgMaterialUtilities::computeMaterialParameters(iMaterialAssignments, iBaseMaterialStructures, oStresses, strains, physicalParameters);
      }    
      writeFiles(iMatStructureSize[0], iMaterialAssignments, iBaseMaterialStructures, physicalParameters, oX, "tmp");
      //writeStressDeformationFile(materials, iMatStructureSize, oStresses, oX);
    }
  }
  //writeStressDeformationFile(iMaterialAssignments, iMatStructureSize, oStresses, oX);
  return 0;
}

void SamplesGeneratorImpl::writeStressDeformationFile(const std::vector<std::vector<int> > &iMaterialAssignments, int iMatStructureSize[3],
                                                      const std::vector<std::vector<float> > iStresses[2], const std::vector<std::vector<std::vector<float> > > iX[2], std::string iPostfix)
{
  int c = m_blockRep;
  int n[3] = {c*iMatStructureSize[0]*m_nbSubdivisions, c*iMatStructureSize[1]*m_nbSubdivisions, m_dim==2? 0: c*iMatStructureSize[2]*m_nbSubdivisions};
  int level = iMatStructureSize[0];
  int iaxis=0;
  for (iaxis=0; iaxis<2; iaxis++)
  {
    Vector3S forceAxis(0,0,0);
    forceAxis[iaxis] = 1.f;

    std::string extension = ".bin"; // ".txt";
    std::string FileName = m_OutputDirectory + "StressDeformation_" + std::to_string(level) + '_' + std::to_string(m_blockRep) + "_" + (iaxis==0? "x": "y") + iPostfix + extension;
    //writeData(FileName, forceAxis, iMaterialAssignments, stresses[iaxis], x[iaxis], n);
    writeDataBinary(FileName, forceAxis, iMaterialAssignments, iStresses[iaxis], iX[iaxis], n);
  }
}

int SamplesGeneratorImpl::computeMaterialParametersLevel1(std::string & iStepperType, bool iWriteSingleFile, bool iWriteFullDeformation)
{
  // level 0
  std::vector<std::vector<int> > baseMaterialStructures;
  baseMaterialStructures.resize(2);
  baseMaterialStructures[0].push_back(0);
  baseMaterialStructures[1].push_back(1);
  int N[3] = {1,1,1};
  int n[3] = {2*m_nbSubdivisions,2*m_nbSubdivisions,2*m_nbSubdivisions};

  // level 1
  int level = 1;
  int nVoxCoarse = pow(2, m_dim);
  int nmat = 2;
  int ncomb=pow(nmat,nVoxCoarse);
  //std::vector<std::vector<int> > materialAssignments(ncomb);
  std::vector<std::vector<int> > materialAssignments;
  bool mirroredStructure = m_orthotropicOnly||m_cubicOnly;
  int icomb;
  for (icomb=0; icomb<ncomb; icomb++)
  {
    //int2Comb(materialAssignments[icomb], icomb, nmat, nVoxCoarse);
    std::vector<int> matAssignment;
    int2Comb(matAssignment, icomb, nmat, nVoxCoarse);

    std::vector<int> cellMaterials;
    getMaterialAssignment(2, 2, matAssignment, N[0], N[1], baseMaterialStructures, 1, 1, m_nbSubdivisions, cellMaterials);

    if (isStructureManifold(n[0], n[1], cellMaterials, N[0], N[1], mirroredStructure, 0))
    {
      materialAssignments.push_back(matAssignment);
    }
 } 

  std::string FileNameMaterialsBaseStructure = m_OutputDirectory + "Materials_" + std::to_string(level)  + ".txt";
  bool ResOk = writeMaterialCombinations(FileNameMaterialsBaseStructure,  baseMaterialStructures);

  computeMaterialParameters(iStepperType, materialAssignments, baseMaterialStructures, N, level, m_blockRep, iWriteSingleFile, iWriteFullDeformation);

  return 0;
}

int SamplesGeneratorImpl::computeMaterialParameters(std::string & iStepperType, int iLevel, int iNbCombPrevLevel, bool iReadSingleFile, bool iWriteSingleFile, bool iReadFullDeformation, bool iWriteFullDeformation)
{
  assert(iLevel>=2);

  int prevLevel = iLevel-1;
  int blockSize = m_blockRep;

  bool ResOk = true;

  std::vector<cfgScalar> physicalParameters; //YoungModulus, density;
  std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;

  std::string materialFile = m_OutputDirectory + "Materials_" + std::to_string(prevLevel) + ".txt";

  std::string baseFileName = "StressStrain_" + std::to_string(prevLevel) + "_" + std::to_string(blockSize);

  if (iReadSingleFile)
  {
    if (iReadFullDeformation)
    {
      std::string baseFileName = "StressDeformation_" + std::to_string(prevLevel) + "_" + std::to_string(blockSize);
      std::string stressDeformationFileNames[2] = {m_OutputDirectory+ baseFileName + "_x.txt", m_OutputDirectory + baseFileName + "_y.txt"};
      ResOk = cfgMaterialUtilities::computeMaterialParametersFromDeformations(materialFile, stressDeformationFileNames, false, physicalParameters, baseMaterialStructures, materialAssignments);
    }
    else
    {
      std::string stressStrainFileNames[2] = {m_OutputDirectory+ baseFileName + "_x.txt", m_OutputDirectory + baseFileName + "_y.txt"};
      ResOk = cfgMaterialUtilities::computeMaterialParameters(materialFile, stressStrainFileNames, physicalParameters, baseMaterialStructures, materialAssignments);
    }
  }
  else
  {
    std::string subDirectories[2];
    subDirectories[0] =  m_OutputDirectory + "x_level" + std::to_string(prevLevel) + "//";
    subDirectories[1] =  m_OutputDirectory + "y_level" + std::to_string(prevLevel) + "//";
    int NbFiles = iNbCombPrevLevel;
    ResOk = cfgMaterialUtilities::computeMaterialParameters(materialFile, subDirectories, baseFileName, NbFiles, physicalParameters, baseMaterialStructures, materialAssignments);
  }

  int sideSize = pow(2, prevLevel);
  int N[3] = {sideSize,sideSize,sideSize};

  std::vector<std::vector<int> > newBaseMaterialStructures;
  int ibase=0, nbase=(int)materialAssignments.size();
  for (ibase=0; ibase<nbase; ibase++)
  {
    std::vector<int> newBaseMaterial;
    getMaterialAssignment(2, 2, materialAssignments[ibase], N[0]/2, N[1]/2, baseMaterialStructures, 1, 1, m_nbSubdivisions ,newBaseMaterial);
    newBaseMaterialStructures.push_back(newBaseMaterial);
  }

  bool useConvexHull = false;
  int naxis = 2;
  int nparam = naxis + 1;
  std::vector<int> convexHull, boundaryVertices, convexHullFaces;
  std::vector<cfgScalar> parameterPoints;
  parameterPoints.insert(parameterPoints.end(), physicalParameters.begin(), physicalParameters.end());
  if (!useConvexHull)
  {
    std::cout << "Computing delaunay triangulation... " << std::endl;
    computeDelaundayTriangulation(parameterPoints, nparam, convexHull, boundaryVertices, convexHullFaces);
    std::cout << "# params = " << parameterPoints.size() << ", # boundary points = " << boundaryVertices.size() << std::endl;
  }
  else
  {
    computeConvexHull(parameterPoints, nparam, boundaryVertices);
  }
  if (0)
  {
    writeVector2File(boundaryVertices, m_OutputDirectory + "convexHull.m");
    writeVector2File(parameterPoints, m_OutputDirectory + "params1.m");
  }

  std::vector<std::vector<int> > clusters;
  clusterByBoundary(N[0], N[1], newBaseMaterialStructures, 2, clusters);

  std::vector<std::vector<int> > newMaterialAssignments;
  bool mirroredStructure = m_orthotropicOnly||m_cubicOnly;

  int version=5;
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
    newBaseMaterialStructures = getSubVector(newBaseMaterialStructures, boundaryVertices);
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
    getNewCombinations(boundaryVertices, physicalParameters, N, newBaseMaterialStructures, nNewComb, newMaterialAssignments, BaseMaterialStructuresIndices);
    //newBaseMaterialStructures = getSubVector(newBaseMaterialStructures, BaseMaterialStructuresIndices);

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
    int nNewComb = pow(nmat,nVoxCoarse);
    int shift = 1000000;
    nNewComb = shift + 1000000;
    int n[2] = {2*N[0]*m_nbSubdivisions, 2*N[1]*m_nbSubdivisions};
    //newMaterialAssignments.resize(nNewComb);
    int icomb;
    //int NbCellToCheck = 1;
    for (icomb=shift; icomb<nNewComb; icomb++)
    {
      std::vector<int> matAssignment;
      int2Comb(matAssignment, icomb, nmat, nVoxCoarse);

      /*if (icomb>nmat*nmat*nmat)
      {
        NbCellToCheck = 4;
      }
      else if (icomb>nmat*nmat)
      {
        NbCellToCheck = 3;
      }
      else if (icomb>nmat)
      {
        NbCellToCheck = 2;
      }*/ 

      std::vector<int> cellMaterials;
      getMaterialAssignment(2, 2, matAssignment, N[0], N[1], newBaseMaterialStructures, m_blockRep, m_blockRep, m_nbSubdivisions, cellMaterials);

      if (isStructureManifold(n[0], n[1], cellMaterials, N[0], N[1], mirroredStructure, 4))
      {
        newMaterialAssignments.push_back(matAssignment);
      }
      //int2Comb(newMaterialAssignments[icomb], icomb, nmat, nVoxCoarse);
    } 
    std::cout << "Nb structure "<< newMaterialAssignments.size() << std::endl; 
  }
  else if (version==3)
  {
    int nNewComb = 1000;
    std::vector<int> BaseMaterialStructuresIndices;
    getNewCombinationsV2(boundaryVertices, physicalParameters, N, nNewComb, newMaterialAssignments, BaseMaterialStructuresIndices);
  }
  else if (version==4)
  {
    //int nNewComb = 10000;
    int nNewComb = 100;
    //getNewCombinationsV3(convexHullFaces, physicalParameters, N, nNewComb, newMaterialAssignments, BaseMaterialStructuresIndices);
    getNewCombinationsV4(convexHullFaces, physicalParameters, N, nNewComb, newMaterialAssignments, newBaseMaterialStructures);
  }
  else if (version==5)
  {
    int nNewComb = 1000000;
    getNewCombinationsV5(convexHullFaces, physicalParameters, N, nNewComb, newMaterialAssignments, newBaseMaterialStructures);
  }
  std::string FileNameMaterials = m_OutputDirectory + "Materials_" + std::to_string(iLevel)  + ".txt";
  writeMaterialCombinations(FileNameMaterials,  newBaseMaterialStructures);

  computeMaterialParameters(iStepperType, newMaterialAssignments, newBaseMaterialStructures, N, iLevel, m_blockRep, iWriteSingleFile, iWriteFullDeformation);

  return 0;
}

void SamplesGeneratorImpl::growStructure(int N[2], const std::vector<std::vector<int> > &materialAssignments, std::vector<std::vector<int> > &oNewMaterialAssignments)
{
  int n[2] = {N[0]+1, N[1]+1};

  std::set<std::vector<int> > newMatAssignments;
  oNewMaterialAssignments.clear();
  int istruct, nstruct=(int)materialAssignments.size();
  for (istruct=0; istruct<nstruct; istruct++)
  {
    const std::vector<int> & matAssignement = materialAssignments[istruct];

    int i;
    for (i=0; i<N[0]; i++)
    {
      std::vector<int> duplicatedLayerX;
      getLayer(N[0], N[1], matAssignement, i, 0, duplicatedLayerX);
      std::vector<int> newMatAssignement1;
      insertLayer(N[0], N[1], matAssignement, duplicatedLayerX, i, 0, newMatAssignement1);

      int j;
      for (j=0; j<N[1]; j++)
      {
        std::vector<int> duplicatedLayerY;
        getLayer(N[0]+1, N[1], newMatAssignement1, j, 1, duplicatedLayerY);
        std::vector<int> newMatAssignement2;
        insertLayer(N[0]+1, N[1], newMatAssignement1, duplicatedLayerY, j, 1, newMatAssignement2);

        newMatAssignments.insert(newMatAssignement2);
      }
    }
  }
  oNewMaterialAssignments = toStdVector(newMatAssignments);
}

void SamplesGeneratorImpl::growStructureDoubleSize(int N[2], const std::vector<std::vector<int> > &materialAssignments, std::vector<std::vector<int> > &oNewMaterialAssignments, int iRowToDuplicate, int iColToDuplicate)
{
  int n[2] = {2*N[0], 2*N[1]};

  std::set<std::vector<int> > newMatAssignments;
  oNewMaterialAssignments.clear();

  int istruct, nstruct=(int)materialAssignments.size();
  for (istruct=0; istruct<nstruct; istruct++)
  {
    const std::vector<int> & matAssignement = materialAssignments[istruct];
    growStructureDoubleSize(N, matAssignement, newMatAssignments, iRowToDuplicate, iColToDuplicate);
  }
  oNewMaterialAssignments = toStdVector(newMatAssignments); 
}



void SamplesGeneratorImpl::growStructureDoubleSize(int N[2], const std::vector<int> &matAssignment, std::set<std::vector<int> > &ioNewMaterialAssignments, int iRowToDuplicate, int iColToDuplicate)
{
  int n[3] = {2*N[0], 2*N[1], 2*N[2]};

  if (matAssignment.size()==1)
  {
    std::vector<int> newMatAssignement(4, matAssignment[0]);
    ioNewMaterialAssignments.insert(newMatAssignement);
  }
  else
  {
    std::vector<std::vector<int> > layersX(N[0]);
    for (int i=0; i<N[0]; i++)
    {
      getLayer(N[0], N[1], matAssignment, i, 0, layersX[i]);
    }

    std::vector<int> newMatAssignement1(n[0]*N[1]);
    int imin=0, imax=3;
    if (iRowToDuplicate>=0)
    {
      imin = iRowToDuplicate;
      imax = iRowToDuplicate+1;
    }
    for (int i=imin; i<imax; i++)
    {
      int scale1 = i+1;
      int scale2 = 3-i;

      int ncol = N[0]/2;
      int nrow = N[1];
      int icol = 0;
      for (icol=0; icol<ncol; icol++)
      {
        for (int irow=0; irow<nrow; irow++)
        {
          int mat = layersX[icol][irow];
          for (int irep=0; irep<scale1; irep++)
          {
            int x = scale1*icol + irep;
            int y = irow;
            int indMat = getGridToVectorIndex(x, y, n[0], N[1]);
            newMatAssignement1[indMat] = mat;
          }
        }
      }
      int shift = scale1*ncol;
      for (icol=0; icol<ncol; icol++)
      {
        for (int irow=0; irow<nrow; irow++)
        {
          int mat = layersX[ncol+icol][irow];
          for (int irep=0; irep<scale2; irep++)
          {
            int x = shift + scale2*icol + irep;
            int y = irow;
            int indMat = getGridToVectorIndex(x, y, n[0], N[1]);
            newMatAssignement1[indMat] = mat;
          }
        }
      }
      std::vector<std::vector<int> > layersY(N[1]);
      for (int j=0; j<N[1]; j++)
      {
        getLayer(n[0], N[1], newMatAssignement1, j, 1, layersY[j]);
      }

      int jmin=0, jmax=3;
      if (iColToDuplicate>=0)
      {
        jmin = iColToDuplicate;
        jmax = iColToDuplicate+1;
      }
      for (int j=jmin; j<jmax; j++)
      {
        std::vector<int> newMatAssignement2(n[0]*n[1]);

        int scale1 = j+1;
        int scale2 = 3-j;

        int ncol = n[0];
        int nrow = N[1]/2;
        int irow = 0;
        for (irow=0; irow<nrow; irow++)
        {
          for (int icol=0; icol<ncol; icol++)
          {
            int mat = layersY[irow][icol];
            for (int irep=0; irep<scale1; irep++)
            {
              int y = scale1*irow + irep;
              int x = icol;
              int indMat = getGridToVectorIndex(x, y, n[0], n[1]);
              newMatAssignement2[indMat] = mat;
            }
          }
        }
        int shift = scale1*nrow;
        for (irow=0; irow<nrow; irow++)
        {
          for (int icol=0; icol<ncol; icol++)
          {
            int mat = layersY[irow+nrow][icol];
            for (int irep=0; irep<scale2; irep++)
            {
              int y = shift + scale2*irow + irep;
              int x = icol;
              int indMat = getGridToVectorIndex(x, y, n[0], n[1]);
              newMatAssignement2[indMat] = mat;
            }
          }
        }
        ioNewMaterialAssignments.insert(newMatAssignement2);
      }
    }
  }
}

void SamplesGeneratorImpl::growStructureDoubleSize(int N[3], const std::vector<int> &matAssignment, std::set<std::vector<int> > &ioNewMaterialAssignments, int iLayerXToDuplicate, int iLayerYToDuplicate, int iLayerZToDuplicate)
{
  int n[3] = {2*N[0], 2*N[1], 2*N[2]};

  if (matAssignment.size()==1)
  {
    std::vector<int> newMatAssignement(8, matAssignment[0]);
    ioNewMaterialAssignments.insert(newMatAssignement);
  }
  else
  {
    std::vector<std::vector<int> > layersX(N[0]);
    for (int i=0; i<N[0]; i++)
    {
      getLayer(N[0], N[1], N[2], matAssignment, i, 0, layersX[i]);
    }

    std::vector<int> newMatAssignement1(n[0]*N[1]*N[2]);
    int imin=0, imax=3;
    if (iLayerXToDuplicate>=0)
    {
      imin = iLayerXToDuplicate;
      imax = iLayerXToDuplicate+1;
    }
    for (int i=imin; i<imax; i++)
    {
      int scale1 = i+1;
      int scale2 = 3-i;

      int nx = N[0]/2;
      int ny = N[1];
      int nz = N[2];
      for (int ix=0; ix<nx; ix++)
      {
        for (int iy=0; iy<ny; iy++)
        {
          for (int iz=0; iz<nz; iz++)
          {
            int mat = layersX[ix][nz*iy+iz];
            for (int irep=0; irep<scale1; irep++)
            {
              int x = scale1*ix + irep;
              int y = iy;
              int z = iz;
              int indMat = getGridToVectorIndex(x, y, z, n[0], N[1], N[2]);
              newMatAssignement1[indMat] = mat;
            }
          }
        }
      }
      int shift = scale1*nx;
      for (int ix=0; ix<nx; ix++)
      {
        for (int iy=0; iy<ny; iy++)
        {
          for (int iz=0; iz<nz; iz++)
          {
            int mat = layersX[nx+ix][nz*iy+iz];
            for (int irep=0; irep<scale2; irep++)
            {
              int x = shift + scale2*ix + irep;
              int y = iy;
              int z = iz;
              int indMat = getGridToVectorIndex(x, y, z, n[0], N[1], N[2]);
              newMatAssignement1[indMat] = mat;
            }
          }
        }
      }
      std::vector<std::vector<int> > layersY(N[1]);
      for (int j=0; j<N[1]; j++)
      {
        getLayer(n[0], N[1], N[2], newMatAssignement1, j, 1, layersY[j]);
      }

      int jmin=0, jmax=3;
      if (iLayerYToDuplicate>=0)
      {
        jmin = iLayerYToDuplicate;
        jmax = iLayerYToDuplicate+1;
      }
      for (int j=jmin; j<jmax; j++)
      {
        std::vector<int> newMatAssignement2(n[0]*n[1]*N[2]);

        int scale1 = j+1;
        int scale2 = 3-j;

        int nx = n[0];
        int ny = N[1]/2;
        for (int iy=0; iy<ny; iy++)
        {
          for (int ix=0; ix<nx; ix++)
          {
            for (int iz=0; iz<nz; iz++)
            {
              int mat = layersY[iy][nz*ix+iz];
              for (int irep=0; irep<scale1; irep++)
              {
                int y = scale1*iy + irep;
                int x = ix;
                int z = iz;
                int indMat = getGridToVectorIndex(x, y, z, n[0], n[1], N[2]);
                newMatAssignement2[indMat] = mat;
              }
            }
          }
        }
        int shift = scale1*ny;
        for (int iy=0; iy<ny; iy++)
        {
          for (int ix=0; ix<nx; ix++)
          {
            for (int iz=0; iz<nz; iz++)
            {
              int mat = layersY[iy+ny][nz*ix+iz];
              for (int irep=0; irep<scale2; irep++)
              {
                int y = shift + scale2*iy + irep;
                int x = ix;
                int z = iz;
                int indMat = getGridToVectorIndex(x, y, z, n[0], n[1], N[2]);
                newMatAssignement2[indMat] = mat;
              }
            }
          }
        }

        std::vector<std::vector<int> > layersZ(N[2]);
        for (int k=0; k<N[2]; k++)
        {
          getLayer(n[0], n[1], N[2], newMatAssignement2, k, 2, layersZ[k]);
        }

        int kmin=0, kmax=3;
        if (iLayerZToDuplicate>=0)
        {
          kmin = iLayerZToDuplicate;
          kmax = iLayerZToDuplicate+1;
        }
        for (int k=kmin; k<kmax; k++)
        {
          std::vector<int> newMatAssignement3(n[0]*n[1]*n[2]);

          int scale1 = k+1;
          int scale2 = 3-k;

          int nx = n[0];
          int ny = n[1];
          int nz = N[2]/2;
          for (int iz=0; iz<nz; iz++)
          {
            for (int ix=0; ix<nx; ix++)
            {
              for (int iy=0; iy<ny; iy++)
              {
                int mat = layersZ[iz][ny*ix+iy];
                for (int irep=0; irep<scale1; irep++)
                {
                  int x = ix;
                  int y = iy;
                  int z = scale1*iz + irep;
                  int indMat = getGridToVectorIndex(x, y, z, n[0], n[1], n[2]);
                  newMatAssignement3[indMat] = mat;
                }
              }
            }
          }
          int shift = scale1*nz;
          for (int iz=0; iz<nz; iz++)
          {
            for (int ix=0; ix<nx; ix++)
            {
              for (int iy=0; iy<ny; iy++)
              {
                int mat = layersZ[iz+nz][ny*ix+iy];
                for (int irep=0; irep<scale2; irep++)
                {
                  int x = ix;
                  int y = iy;
                  int z = shift + scale2*iz + irep;
                  int indMat = getGridToVectorIndex(x, y, z, n[0], n[1], n[2]);
                  newMatAssignement3[indMat] = mat;
                }
              }
            }
          }
          ioNewMaterialAssignments.insert(newMatAssignement3);
        }
      }
    }
  }
}

int SamplesGeneratorImpl::computeMaterialParametersIncremental(std::string & iStepperType, int iLevel, std::vector<std::vector<int> > &oNewMaterialAssignments)
{
  bool ResOk = true;

  std::vector<std::vector<int> > newMaterialAssignments;
  std::vector<std::vector<int> > newBaseMaterialStructures;

  int n[3] = {iLevel, iLevel, iLevel};
  int N[3] = {1,1,1};

  if (iLevel<2)
  {
    std::vector<int> mat1(1,0), mat2(1,1);
    newBaseMaterialStructures.resize(2);
    newBaseMaterialStructures[0] = mat1;
    newBaseMaterialStructures[1] = mat2;

    newMaterialAssignments = newBaseMaterialStructures;
  }
  else
  {
    int readBinary = true;
    bool doubleSize = true;

    int prevLevel = doubleSize? iLevel/2: iLevel-1;
    int blockSize = m_blockRep;

    std::vector<cfgScalar> physicalParameters; //YoungModulus, density;
    std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;

    ResOk = readFiles(prevLevel, materialAssignments, baseMaterialStructures, physicalParameters);
  
    int n_prev[3] = {prevLevel, prevLevel, prevLevel};  
    newBaseMaterialStructures = baseMaterialStructures;

    if (m_cubicOnly)
    {
      int nprev[3] = {prevLevel, prevLevel, prevLevel};
      int nVoxCoarse = getDiagonalStructureNumberOfElements(prevLevel, m_dim);
      int nmat = 2;
      int ncomb=pow(nmat,nVoxCoarse);
      std::cout << "** Nb comb = " << ncomb << std::endl;
      for (int icomb=0; icomb<ncomb; icomb++)
      {
        std::vector<int> matAssignment;
        int2Comb(matAssignment, icomb, nmat, nVoxCoarse);
        if (m_dim==2)
        {
          std::vector<int> matAssignmentMirroredAlongDiag;
          mirrorStructureAlongDiagonal(nprev[0], nprev[1], matAssignment, matAssignmentMirroredAlongDiag);
          if (isStructureManifold(nprev[0], nprev[1], matAssignmentMirroredAlongDiag, N[0], N[1], true, 0))
          {
            std::vector<int> matAssignmentMirrored;
            mirrorStructure(nprev[0], nprev[1], matAssignmentMirroredAlongDiag, matAssignmentMirrored);
            newMaterialAssignments.push_back(matAssignmentMirrored);
          }
        }
        else
        {
          std::vector<int> matAssignmentMirroredAlongDiag;
          mirrorStructureAlongDiagonal(nprev[0], nprev[1], nprev[2], matAssignment, matAssignmentMirroredAlongDiag);
          if (isStructureManifold(nprev[0], nprev[1], nprev[2], matAssignmentMirroredAlongDiag, N[0], N[1], N[2], true))
          {
            std::vector<int> matAssignmentMirrored;
            mirrorStructure(nprev[0], nprev[1], nprev[2], matAssignmentMirroredAlongDiag, matAssignmentMirrored);
            newMaterialAssignments.push_back(matAssignmentMirrored);
          }
        }
      } 
    }
    else if (m_orthotropicOnly)
    {
      int nprev[3] = {prevLevel, prevLevel, prevLevel};
      int nVoxCoarse = pow(prevLevel, m_dim);
      int nmat = 2;
      int ncomb=pow(nmat,nVoxCoarse);
      std::cout << "** Nb comb = " << ncomb << std::endl;
      for (int icomb=0; icomb<ncomb; icomb++)
      {
        std::vector<int> matAssignment;
        int2Comb(matAssignment, icomb, nmat, nVoxCoarse);
        if (m_dim==2)
        {
          if (isStructureManifold(nprev[0], nprev[1], matAssignment, N[0], N[1], true, 0))
          {
            std::vector<int> matAssignmentMirrored;
            mirrorStructure(nprev[0], nprev[1], matAssignment, matAssignmentMirrored);
            newMaterialAssignments.push_back(matAssignmentMirrored);
          }
        }
        else
        {
          if (isStructureManifold(nprev[0], nprev[1], nprev[2], matAssignment, N[0], N[1], N[2], true))
          {
            std::vector<int> matAssignmentMirrored;
            mirrorStructure(nprev[0], nprev[1], nprev[2], matAssignment, matAssignmentMirrored);
            newMaterialAssignments.push_back(matAssignmentMirrored);
          }
        }
      } 
    }
    else
    {
      if (iLevel==2)
      {
        int nVoxCoarse = pow(2, m_dim);
        int nmat = 2;
        int ncomb=pow(nmat,nVoxCoarse);
        int icomb;
        for (icomb=0; icomb<ncomb; icomb++)
        {
          std::vector<int> matAssignment;
          int2Comb(matAssignment, icomb, nmat, nVoxCoarse);
          if (m_dim==2)
          {
            if (isStructureManifold(n[0], n[1], matAssignment, N[0], N[1], 0, false))
            {
              newMaterialAssignments.push_back(matAssignment);
            }
          }
          else
          {
            if (isStructureManifold(n[0], n[1], n[2], matAssignment, N[0], N[1], N[2], false))
            {
              newMaterialAssignments.push_back(matAssignment);
            }
          }
        } 
      }
      else
      {
        if (doubleSize)
        {
          growStructureDoubleSize(n_prev,  materialAssignments, newMaterialAssignments);
        }
        else
        {
          growStructure(n_prev, materialAssignments, newMaterialAssignments);
        }
      }
    }
  }
  std::cout << "** Nb manifold comb = " << newMaterialAssignments.size() << std::endl;
  oNewMaterialAssignments = newMaterialAssignments;

  return ResOk;
}

int SamplesGeneratorImpl::computeVariations(std::string & iStepperType, int iLevel)
{
  bool ResOk = true;

  std::vector<std::vector<int> > newBaseMaterialStructures;

  int n[3] = {iLevel, iLevel, iLevel};
  int N[3] = {1,1,1};

  std::vector<std::vector<std::vector<int> > > materialAssignments, baseMaterialStructures;
  std::vector<std::vector<float> > physicalParameters;
  int ilevel=0, nlevel=iLevel;
  materialAssignments.resize(nlevel+1);
  baseMaterialStructures.resize(nlevel+1);
  physicalParameters.resize(nlevel+1);
  for (ilevel=1; ilevel<=nlevel; ilevel++)
  {
    ResOk = readFiles(ilevel, materialAssignments[ilevel], baseMaterialStructures[ilevel], physicalParameters[ilevel]);
  }

  std::vector<std::multimap<float, int> > dist2Mat(nlevel);
  dist2Mat[1].insert(std::pair<float,int>(0,0));
  dist2Mat[1].insert(std::pair<float,int>(0,1));

  for (ilevel=2; ilevel<nlevel; ilevel++)
  {
    std::vector<int> convexHull, boundaryVertices, boundaryFaces;
    computeDelaundayTriangulation(physicalParameters[ilevel], 3, convexHull, boundaryVertices, boundaryFaces);

    std::vector<cfgScalar> & points = physicalParameters[ilevel];
    std::vector<cfgScalar> sampledBoundary;
    sampleMesh(points, boundaryFaces, 1000, sampledBoundary);

    std::vector<cfgScalar> boundaryPoints = getSubVector(points, 3, boundaryVertices);
    sampledBoundary.insert(sampledBoundary.end(), boundaryPoints.begin(), boundaryPoints.end());

    std::vector<cfgScalar> distToBoundary;
    computeDistancesToPointCloud(points, sampledBoundary, distToBoundary);
    cfgScalar maxDist = *std::max_element(distToBoundary.begin(), distToBoundary.end()); 
    distToBoundary = cfgUtil::mult<cfgScalar>(distToBoundary, 1./maxDist);

    int ipoint=0, npoint=(int)points.size()/3;
    for (ipoint=0; ipoint<npoint; ipoint++)
    {
      dist2Mat[ilevel].insert(std::pair<float,int>(distToBoundary[ipoint],ipoint));
    }
  }
  bool mirroredStructure = m_orthotropicOnly||m_cubicOnly;
  std::vector<std::vector<int> > newMaterialAssignments;
  int level = nlevel;
  int istruct=0, nstruct=(int)materialAssignments[level].size();
  std::cout << "Nb structures = " << nstruct << std::endl;
  for (istruct=0; istruct<nstruct; istruct++)
  {
    int structIndex = istruct;
    //int structIndex = 295;
    //int structIndex = 242; 
    //int structIndex = 222;
  
    int nsample = 10;
    float eps = 1.e-4;
    std::vector<int> initMaterialAsssignment = materialAssignments[level][structIndex];
    int isample=0;
    while (isample<nsample)
    {
      std::vector<int> materialAsssignment = initMaterialAsssignment;
      int nvoxel = materialAsssignment.size();

      int voxelIndex = rand() % nvoxel;
      int subLevel = 8;
      int nsub[3] = {subLevel, subLevel, subLevel};
      int matIndex = -1;
      if (subLevel > 1)
      {
        float r = (float)rand()/float(RAND_MAX);
        matIndex = getRandomMatWithLinearProba(dist2Mat[subLevel], eps, r);
      }
      else
      {
        matIndex = rand() % materialAssignments[subLevel].size();
      }
      if (matIndex >= 0)
      {
        std::vector<int> &subMaterialStructure = materialAssignments[subLevel][matIndex];
        updateMaterialSubstructure(materialAsssignment, n[0], n[1], subMaterialStructure, nsub[0], nsub[1], voxelIndex);
        if (materialAsssignment != initMaterialAsssignment && isStructureManifold(n[0], n[1], materialAsssignment, mirroredStructure))
        {
          newMaterialAssignments.push_back(materialAsssignment);
          isample++;
        } 
      }
    }
  }
  if (0)
  {
    std::vector<int> materialAsssignmentCol1, materialAsssignmentCol2, materialAsssignmentCol3;
    materialAsssignmentCol1.resize(10, 0);
    materialAsssignmentCol1[0]=1;
    materialAsssignmentCol1[5]=1;

    materialAsssignmentCol2.resize(10, 1);

    materialAsssignmentCol3.resize(10, 0);
    materialAsssignmentCol3[3]=1;
    materialAsssignmentCol3[8]=1;

     std::vector<int> newMaterialAssignment;
     newMaterialAssignment.insert(newMaterialAssignment.end(), materialAsssignmentCol1.begin(), materialAsssignmentCol1.end());
     newMaterialAssignment.insert(newMaterialAssignment.end(), materialAsssignmentCol2.begin(), materialAsssignmentCol2.end());
     newMaterialAssignment.insert(newMaterialAssignment.end(), materialAsssignmentCol2.begin(), materialAsssignmentCol2.end());
     newMaterialAssignment.insert(newMaterialAssignment.end(), materialAsssignmentCol3.begin(), materialAsssignmentCol3.end());
     newMaterialAssignment.insert(newMaterialAssignment.end(), materialAsssignmentCol2.begin(), materialAsssignmentCol2.end());
     newMaterialAssignment.insert(newMaterialAssignment.end(), materialAsssignmentCol1.begin(), materialAsssignmentCol1.end());
     newMaterialAssignment.insert(newMaterialAssignment.end(), materialAsssignmentCol2.begin(), materialAsssignmentCol2.end());
     newMaterialAssignment.insert(newMaterialAssignment.end(), materialAsssignmentCol3.begin(), materialAsssignmentCol3.end());
     newMaterialAssignment.insert(newMaterialAssignment.end(), materialAsssignmentCol2.begin(), materialAsssignmentCol2.end());
     newMaterialAssignment.insert(newMaterialAssignment.end(), materialAsssignmentCol1.begin(), materialAsssignmentCol1.end());
     newMaterialAssignments.push_back(newMaterialAssignment);
  }
  if (0)
  {
    std::vector<std::vector<int> >materialAsssignmentPerCol(10);
    for (int icol=0; icol<10; icol++)
    {
      materialAsssignmentPerCol[icol].resize(10, 1);
    }
    materialAsssignmentPerCol[0][0] = 0;
    materialAsssignmentPerCol[1][0] = 0;
    materialAsssignmentPerCol[2][0] = 0;
    materialAsssignmentPerCol[3][0] = 0;
    materialAsssignmentPerCol[4][0] = 0;

    materialAsssignmentPerCol[5][5] = 0;
    materialAsssignmentPerCol[6][5] = 0;
    materialAsssignmentPerCol[7][5] = 0;
    materialAsssignmentPerCol[8][5] = 0;
    materialAsssignmentPerCol[9][5] = 0;

    materialAsssignmentPerCol[2][3] = 0;
    materialAsssignmentPerCol[2][4] = 0;
    materialAsssignmentPerCol[2][5] = 0;
    materialAsssignmentPerCol[2][6] = 0;
    materialAsssignmentPerCol[2][7] = 0;

    materialAsssignmentPerCol[7][0] = 0;
    materialAsssignmentPerCol[7][1] = 0;
    materialAsssignmentPerCol[7][2] = 0;
    materialAsssignmentPerCol[7][8] = 0;
    materialAsssignmentPerCol[7][9] = 0;

    std::vector<int> newMaterialAssignment;
    for (int icol=0; icol<10; icol++)
    {
      newMaterialAssignment.insert(newMaterialAssignment.end(), materialAsssignmentPerCol[icol].begin(), materialAsssignmentPerCol[icol].end());
    }
    newMaterialAssignments.push_back(newMaterialAssignment);
  }
  computeMaterialParameters(iStepperType, newMaterialAssignments, n, baseMaterialStructures[level], N, level, m_blockRep, true, true, "variations");
  return ResOk;
}

int SamplesGeneratorImpl::computeVariationsV2(int iLevel, int iSubLevel, 
                                              const std::vector<std::vector<std::vector<int> > > &iMaterialAssignmentsPreviousLevels,  
                                              const std::vector<std::vector<std::vector<int> > > &iBaseMaterialStructuresPreviousLevels, 
                                              const std::vector<std::vector<float> > &iPhysicalParametersPreviousLevels,
                                              std::vector<std::vector<int> > &oNewMaterialAssignments)
{
  oNewMaterialAssignments.clear();

  bool ResOk = true;

  std::vector<std::vector<int> > newBaseMaterialStructures;

  int n[3] = {iLevel, iLevel, iLevel};
  int N[3] = {1,1,1};

  int nbClusters = 200;
  std::vector<cfgScalar> lengths(3, 1);
  int nbIter = 10;

  const std::vector<std::vector<std::vector<int> > > & materialAssignments = iMaterialAssignmentsPreviousLevels;
  const std::vector<std::vector<std::vector<int> > > & baseMaterialStructures = iBaseMaterialStructuresPreviousLevels;
  const std::vector<std::vector<float> > & physicalParameters = iPhysicalParametersPreviousLevels;

  int ilevel=0, nlevel=iLevel;
  std::vector<std::vector<std::vector<int> > > clusters;
  clusters.resize(nlevel+1);
  for (ilevel=1; ilevel<=nlevel; ilevel*=2)
  {
    std::vector<float> rescaledParameters = physicalParameters[ilevel];
    rescaleData(rescaledParameters, 3, lengths);
    std::vector<int> pointIndices; 
    getKMeans(nbIter, nbClusters, rescaledParameters, 3, clusters[ilevel], &pointIndices);
    if (1)
    {
      std::vector<std::vector<std::vector<float> > > clusted_x[2];
      std::vector<float> clusteredPhysicalParameters = getSubVector(physicalParameters[ilevel], 3, pointIndices);
      std::vector<std::vector<int> > clusteredMaterialAssignments = getSubVector(materialAssignments[ilevel], pointIndices);

      writeFiles(ilevel, clusteredMaterialAssignments, baseMaterialStructures[ilevel], clusteredPhysicalParameters, clusted_x, "clustered");
    }
  }

  int level = nlevel;
  bool mirroredStructure = m_orthotropicOnly||m_cubicOnly;
  int icluster=0, ncluster=(int)clusters[level].size();
  std::cout << "Nb clusters = " << ncluster << std::endl;
  for (icluster=0; icluster<ncluster; icluster++)
  {
    int itrial=0, ntrial=5;
    for (itrial=0; itrial<ntrial; itrial++)
    {
      int indPointInCluster = rand() % clusters[level][icluster].size();
      int structIndex = clusters[level][icluster][indPointInCluster];

      int nsample = 5;
      float eps = 1.e-4;
      std::vector<int> initMaterialAsssignment = materialAssignments[level][structIndex];
      int isample=0;
      while (isample<nsample)
      {
        std::vector<int> materialAsssignment = initMaterialAsssignment;
        int nvoxel = materialAsssignment.size();

        int voxelIndex = rand() % nvoxel;
        int subLevel = iSubLevel;
        int nsub[3] = {subLevel, subLevel, subLevel};

        int indCluster = rand() % clusters[subLevel].size();
        int indMat = rand() % clusters[subLevel][indCluster].size();
        int matIndex = clusters[subLevel][indCluster][indMat];

        if (matIndex >= 0)
        {
          const std::vector<int> &subMaterialStructure = materialAssignments[subLevel][matIndex];
          updateMaterialSubstructure(materialAsssignment, n[0], n[1], subMaterialStructure, nsub[0], nsub[1], voxelIndex);
          if (materialAsssignment != initMaterialAsssignment && isStructureManifold(n[0], n[1], materialAsssignment, mirroredStructure))
          {
            oNewMaterialAssignments.push_back(materialAsssignment);
            isample++;
          } 
        }
      }
    }
  }
  return ResOk;
}

int SamplesGeneratorImpl::computeVariationsV3(int iLevel, int iSubLevel, 
                                              const std::vector<std::vector<std::vector<int> > > &iMaterialAssignmentsPreviousLevels,  
                                              const std::vector<std::vector<std::vector<int> > > &iBaseMaterialStructuresPreviousLevels, 
                                              const std::vector<std::vector<float> > &iPhysicalParametersPreviousLevels,
                                              std::vector<std::vector<int> > &oNewMaterialAssignments)
{
  oNewMaterialAssignments.clear();

  bool ResOk = true;

  std::vector<std::vector<int> > newBaseMaterialStructures;

  int n[3] = {iLevel, iLevel, iLevel};
  int N[3] = {1,1,1};

  int nbClusters = 200;
  std::vector<cfgScalar> lengths(3, 1);
  int nbIter = 10;

  const std::vector<std::vector<std::vector<int> > > & materialAssignments = iMaterialAssignmentsPreviousLevels;
  const std::vector<std::vector<std::vector<int> > > & baseMaterialStructures = iBaseMaterialStructuresPreviousLevels;
  const std::vector<std::vector<float> > & physicalParameters = iPhysicalParametersPreviousLevels;

  int ilevel=0, nlevel=iLevel;
  std::vector<std::vector<std::vector<int> > > clusters;
  clusters.resize(nlevel+1);
  for (ilevel=1; ilevel<nlevel; ilevel*=2)
  {
    std::vector<float> rescaledParameters = physicalParameters[ilevel];
    rescaleData(rescaledParameters, 3, lengths);
    std::vector<int> pointIndices; 
    getKMeans(nbIter, nbClusters, rescaledParameters, 3, clusters[ilevel], &pointIndices);
    if (1)
    {
      std::vector<std::vector<std::vector<float> > > clusted_x[2];
      std::vector<float> clusteredPhysicalParameters = getSubVector(physicalParameters[ilevel], 3, pointIndices);
      std::vector<std::vector<int> > clusteredMaterialAssignments = getSubVector(materialAssignments[ilevel], pointIndices);

      writeFiles(ilevel, clusteredMaterialAssignments, baseMaterialStructures[ilevel], clusteredPhysicalParameters, clusted_x, "clustered");
    }
  }

  int level = nlevel/2;
  bool mirroredStructure = m_orthotropicOnly||m_cubicOnly;
  int icluster=0, ncluster=(int)clusters[level].size();
  std::cout << "Nb clusters = " << ncluster << std::endl;
  for (icluster=0; icluster<ncluster; icluster++)
  {
    int itrial=0, ntrial=5;
    for (itrial=0; itrial<ntrial; itrial++)
    {
      int indPointInCluster = rand() % clusters[level][icluster].size();
      int structIndex = clusters[level][icluster][indPointInCluster];

      std::vector<int> materialAsssignmentPrevLevel = materialAssignments[level][structIndex];
      std::set<std::vector<int> > initMaterialAsssignments;
      int prevLevel = ilevel/2;
      int n_prev[3] = {prevLevel, prevLevel, prevLevel};  
      int rowToDuplicate = rand()%3;
      int colToDuplicate = rand()%3;
      growStructureDoubleSize(n_prev,  materialAsssignmentPrevLevel, initMaterialAsssignments, rowToDuplicate, colToDuplicate);
      std::vector<int> initMaterialAsssignment = *initMaterialAsssignments.begin();

      int nsample = 5;
      float eps = 1.e-4;
      int isample=0;
      while (isample<nsample)
      {
        std::vector<int> materialAsssignment = initMaterialAsssignment;
        int nvoxel = materialAsssignment.size();

        int voxelIndex = rand() % nvoxel;
        int subLevel = iSubLevel;
        int nsub[3] = {subLevel, subLevel, subLevel};

        int indCluster = rand() % clusters[subLevel].size();
        int indMat = rand() % clusters[subLevel][indCluster].size();
        int matIndex = clusters[subLevel][indCluster][indMat];

        if (matIndex >= 0)
        {
          const std::vector<int> &subMaterialStructure = materialAssignments[subLevel][matIndex];
          updateMaterialSubstructure(materialAsssignment, n[0], n[1], subMaterialStructure, nsub[0], nsub[1], voxelIndex);
          if (materialAsssignment != initMaterialAsssignment && isStructureManifold(n[0], n[1], materialAsssignment, mirroredStructure))
          {
            oNewMaterialAssignments.push_back(materialAsssignment);
            isample++;
          } 
        }
      }
    }
  }
  return ResOk;
}

int SamplesGeneratorImpl::getNbParameters()
{
  int nparam = 3;
  if (m_cubicOnly)
  {
    nparam = 4;  //density, Y, nu
  }
  else if (m_orthotropicOnly)
  {
    if (m_dim==2)
    {
      nparam = 5;  // density, Yx, Yy, Nu_xy, mu
    }
    else if (m_dim==3)
    {
      nparam = 10; // density, Yx, Yy, Yz, Nu_xy, Nu_xz, Nu_yz, mu_x, mu_y, mu_z
    }
  }
  else
  {
    //nparam  = 21; // entries en C tensor
  }
  return nparam;
}

int SamplesGeneratorImpl::getNbFreeCells(int n[3])
{
  int nbCells = n[0]*n[1]*n[2];
  if (m_cubicOnly)
  {
    nbCells = getDiagonalStructureNumberOfElements(n[0]/2, m_dim);
  }
  else if (m_orthotropicOnly)
  {
    if (m_dim==2)
    {
      nbCells/=4;
    }
    else if (m_dim==3)
    {
      nbCells/=8;
    }
  }
  return nbCells;
}

std::vector<int> SamplesGeneratorImpl::getFreeCellsMaterialAssignment(int n[3], const std::vector<int> &iMatAssignments)
{
  std::vector<int> matAssignments;
  if (m_cubicOnly)
  {
    if (m_dim==2)
    {
      getTriangularStructure(n[0], n[1], iMatAssignments, matAssignments);
    }
    else
    {
      getTetrahedralStructure(n[0], n[1], n[2], iMatAssignments, matAssignments);
    }
  }
  else if (m_orthotropicOnly)
  {
    if (m_dim==2)
    {
      getQuarter(n[0], n[1], iMatAssignments, matAssignments);
    }
    else
    {
      getQuarter(n[0], n[1], n[2], iMatAssignments, matAssignments);
    }
  }
  else
  {
    matAssignments = iMatAssignments;
  }
  return matAssignments;
}

std::vector<int>  SamplesGeneratorImpl::getFullMaterialAssignment(int n[3], const std::vector<int> &iFreeCellMatAssignments)
{
  std::vector<int> matAssignments;
  if (m_cubicOnly)
  {
    if (m_dim==2)
    {
      std::vector<int> matAssignmentMirroredAlongDiag;
      mirrorStructureAlongDiagonal(n[0]/2, n[1]/2, iFreeCellMatAssignments, matAssignmentMirroredAlongDiag);
      mirrorStructure(n[0]/2, n[1]/2, matAssignmentMirroredAlongDiag, matAssignments);

      //mirrorStructure(n[0]/2, n[1]/2, iFreeCellMatAssignments, matAssignments);
    }
    else
    {
      std::vector<int> matAssignmentMirroredAlongDiag;
      mirrorStructureAlongDiagonal(n[0]/2, n[1]/2, n[2]/2, iFreeCellMatAssignments, matAssignmentMirroredAlongDiag);
      mirrorStructure(n[0]/2, n[1]/2, n[2]/2, matAssignmentMirroredAlongDiag, matAssignments);

      //mirrorStructure(n[0]/2, n[1]/2, n[2]/2, iFreeCellMatAssignments, matAssignments);
    }
  }
  else if (m_orthotropicOnly)
  {
    if (m_dim==2)
    {
       mirrorStructure(n[0]/2, n[1]/2, iFreeCellMatAssignments, matAssignments);
    }
    else
    {
      mirrorStructure(n[0]/2, n[1]/2, n[2]/2, iFreeCellMatAssignments, matAssignments);
    }
  }
  else
  {
    matAssignments = iFreeCellMatAssignments;
  }
  return matAssignments;
}

void SamplesGeneratorImpl::runContinuousOptimization(int iLevel, 
                                                     const std::vector<MaterialQuad2D> &iBaseMaterials,
                                                     const std::vector<std::vector<int> > &iMaterialAssignments, 
                                                     const std::vector<float> &iParameters,
                                                     const std::vector<float> &iTensors,
                                                     std::vector<std::vector<int> > &oNewMaterialAssignments,
                                                     std::vector<cfgScalar> &oNewPhysicalParameters,
                                                     std::vector<cfgScalar> &oNewTensors)
{
  oNewMaterialAssignments.clear();
  oNewPhysicalParameters.clear();
  oNewTensors.clear();

  int paramdim = getNbParameters();
  DistanceField distanceField(paramdim);
  std::vector<cfgScalar> derivatives;
  std::cout << "Computing distances..." << std::endl;
  std::vector<cfgScalar> distances = distanceField.computeDistances(iParameters, &derivatives);
  cfgScalar coef = 1;
  std::vector<cfgScalar> newPoints = cfgUtil::add<cfgScalar>(iParameters, cfgUtil::mult<cfgScalar>(derivatives, coef));

  int nTargetParticules = 100;
  cfgScalar minRadius = 0.1;
  std::vector<int> newparticules;
  Resampler resampler;
  std::cout << "Resampling boundary..." << std::endl;
  resampler.resampleBoundary(minRadius, paramdim, iParameters, distances, nTargetParticules, newparticules);

  int icycle = 0;
  std::vector<std::vector<int> > baseMaterialStructures(2);
  baseMaterialStructures[0].push_back(0);
  baseMaterialStructures[1].push_back(1);

  std::vector<cfgScalar> initParameters, initTensors;
  initParameters = getSubVector(iParameters, paramdim, newparticules);
  initTensors = getSubVector(iTensors, m_dim==2?6:21, newparticules);
  std::vector<std::vector<int> > &initMaterialAssignments = getSubVector(iMaterialAssignments, newparticules);
  writeFiles(iLevel, initMaterialAssignments, baseMaterialStructures, initParameters, initTensors, "ContinuousOptimizationInit_"+std::to_string(icycle));

  MaterialOptimizer optimizer;
  optimizer.setBaseMaterials(iBaseMaterials);
  if (m_cubicOnly)
  {
    optimizer.setStructureType(MaterialOptimizer::Cubic);
  }
  else if (m_orthotropicOnly)
  {
    optimizer.setStructureType(MaterialOptimizer::Orthotropic);
  }
  std::cout << "Running continuous optimization..." << std::endl;
  int N[2] = {iLevel, iLevel};
  std::vector<std::vector<int> > newMaterialAssignments;
  int istruct=0, nstruct=(int)newparticules.size();
  for (istruct=0; istruct<nstruct; istruct++)
  {
    std::cout << "mat #" << istruct << std::endl;
    int indStruct = newparticules[istruct];
    std::vector<int> matAssignment = iMaterialAssignments[indStruct];
    std::vector<float> targetParameters;
    for (int icoord=0; icoord<paramdim; icoord++)
    {
      targetParameters.push_back(newPoints[paramdim*indStruct+icoord]);
    }
    bool resOk = optimizer.run(N, matAssignment, targetParameters); 
    if (resOk)
    {
      newMaterialAssignments.push_back(matAssignment);
    }
  }

  std::vector<cfgScalar> newParameters, newTensors;
  computeParametersAndTensorValues(N, newMaterialAssignments, newParameters, newTensors);
  writeFiles(iLevel, newMaterialAssignments, baseMaterialStructures, newParameters, newTensors, "ContinousOptimizationResult_"+std::to_string(icycle));
}


void SamplesGeneratorImpl::computeVariationsSMC(int iLevel, 
                                               std::string & iStepperType,
                                               const std::vector<std::vector<std::vector<int> > > &iMaterialAssignmentsPreviousLevels,  
                                               const std::vector<std::vector<std::vector<int> > > &iBaseMaterialStructuresPreviousLevels, 
                                               const std::vector<std::vector<float> > &iPhysicalParametersPreviousLevels,
                                               const std::vector<std::vector<float> > &iTensorsPreviousLevels,
                                               std::vector<std::vector<int> > &oNewMaterialAssignments,
                                               std::vector<cfgScalar> &oNewPhysicalParameters,
                                               std::vector<cfgScalar> &oNewTensors)
{
  oNewMaterialAssignments.clear();

  bool ResOk = true;

  int prevLevel = iLevel/2;
  int n_prev[3] = {prevLevel, prevLevel, (m_dim==3? prevLevel: 1)};  
  int n[3] = {iLevel, iLevel, (m_dim==3? iLevel: 1)};
  int N[3] = {1,1,1};

  const std::vector<std::vector<std::vector<int> > > & baseMaterialStructures = iBaseMaterialStructuresPreviousLevels;

  std::vector<float> physicalParameters = iPhysicalParametersPreviousLevels[prevLevel];
  std::vector<float> tensors = iTensorsPreviousLevels[prevLevel];

  int nTargetParticules = 300;
  int dimparam = getNbParameters();
  int nInitParticules = (int)physicalParameters.size()/dimparam;

  // Sampling
  std::vector<int> particules;
  ScoringFunction scoringFunction(dimparam);
  //cfgScalar radius = scoringFunction.estimateDensityRadius(physicalParameters);
  bool useDistanceField = true;
  scoringFunction.setUseDistanceField(useDistanceField);
  std::vector<cfgScalar> scores = scoringFunction.computeScores(physicalParameters);
  Resampler resampler;
  resampler.resample(scores, nTargetParticules, particules);
 
  int nparticule = (int)particules.size();
  std::vector<std::vector<int> > voxelsToProcess(nparticule);

  int nvoxel = getNbFreeCells(n);
  std::vector<int> initVoxelIndices = genIncrementalSequence(0, nvoxel-1);

  std::vector<std::vector<int> > materialAssignments, allMaterialAssignments;
  std::vector<float> allPhysicalParameters;
  std::vector<float> allTensors;

  bool growStructure = true;
  if (growStructure)
  {
    int istruct=0, nstruct = (int)iMaterialAssignmentsPreviousLevels[prevLevel].size();
    for (istruct=0; istruct<nstruct; istruct++)
    {
      std::vector<int> materialAsssignmentPrevLevel = iMaterialAssignmentsPreviousLevels[prevLevel][istruct];
      std::set<std::vector<int> > initMaterialAsssignments;
      if (m_dim==2)
      {
        growStructureDoubleSize(n_prev,  materialAsssignmentPrevLevel, initMaterialAsssignments, 1, 1);
      }
      else
      {
        growStructureDoubleSize(n_prev,  materialAsssignmentPrevLevel, initMaterialAsssignments, 1, 1, 1);
      }
      std::vector<int> initMaterialAssignment = *initMaterialAsssignments.begin();
      materialAssignments.push_back(initMaterialAssignment);
    }
    allMaterialAssignments = materialAssignments;
    materialAssignments = getSubVector(materialAssignments, 1, particules);
  }
  //allMaterialAssignments = materialAssignments;
  allPhysicalParameters = physicalParameters;
  allTensors = tensors;

  for (int iparticule=0; iparticule<nparticule; iparticule++)
  {
    voxelsToProcess[iparticule] = initVoxelIndices;

    /*int indStructure = particules[iparticule];
    std::vector<int> materialAsssignmentPrevLevel = iMaterialAssignmentsPreviousLevels[prevLevel][indStructure];
    std::set<std::vector<int> > initMaterialAsssignments;
    if (m_dim==2)
    {
      growStructureDoubleSize(n_prev,  materialAsssignmentPrevLevel, initMaterialAsssignments, 1, 1);
    }
    else
    {
      growStructureDoubleSize(n_prev,  materialAsssignmentPrevLevel, initMaterialAsssignments, 1, 1, 1);
    }
    std::vector<int> initMaterialAssignment = *initMaterialAsssignments.begin();
    materialAssignments.push_back(initMaterialAssignment);*/ 
  }
  physicalParameters = getSubVector(physicalParameters, dimparam, particules);
  tensors = getSubVector(tensors, m_dim==2?6:21, particules);

  writeFiles(iLevel, allMaterialAssignments, baseMaterialStructures[prevLevel], allPhysicalParameters, allTensors, "SMC_init");

  std::vector<std::vector<std::vector<float> > > x[3];
  //writeFiles(iLevel, materialAssignments, baseMaterialStructures[prevLevel], physicalParameters, x, "SMC_init");

  bool mirroredStructure = m_orthotropicOnly||m_cubicOnly;

  std::set<int> particulesToProcess = toStdSet(genIncrementalSequence(0, nparticule-1));
  int icycle=0;
  while (particulesToProcess.size())
  {
    std::vector<std::vector<int> > newVoxelsToProcess;
    std::vector<std::vector<int> > newMaterialAssignments;
    std::set<int>::iterator it, it_end=particulesToProcess.end();
    for (it=particulesToProcess.begin(); it!=it_end; it++)
    {
      int indParticule = *it;
 
      std::vector<int> newMaterialAssignment = getFreeCellsMaterialAssignment(n, materialAssignments[indParticule]);

      bool newParticuleGenerated = false;
      while (voxelsToProcess[indParticule].size() && !newParticuleGenerated)
      {
        int randomIndex = rand() % voxelsToProcess[indParticule].size();
        int indVoxel = voxelsToProcess[indParticule][randomIndex];
        voxelsToProcess[indParticule].erase(voxelsToProcess[indParticule].begin()+randomIndex);

        int newMat = rand()%2;

        int oldMat = newMaterialAssignment[indVoxel];
        if (oldMat != newMat)
        {
          newMaterialAssignment[indVoxel] = newMat;
          if (m_cubicOnly)
          {
            std::vector<int> mirroredMaterials;
            if (m_dim==2)
            {
              mirrorStructureAlongDiagonal(n[0]/2, n[1]/2, newMaterialAssignment, mirroredMaterials);
            }
            else
            {
              mirrorStructureAlongDiagonal(n[0]/2, n[1]/2, n[2]/2, newMaterialAssignment, mirroredMaterials);
            }
            if ((m_dim==2 && isStructureManifold(n[0]/2, n[1]/2, mirroredMaterials, mirroredStructure)) || (m_dim==3 && isStructureManifold(n[0]/2, n[1]/2, n[2]/2, mirroredMaterials, N[0], N[1], N[2], mirroredStructure)) )
            {
              newParticuleGenerated = true;
            }
            else
            {
              newMaterialAssignment[indVoxel] = oldMat;
            }
          }
          else
          {
            if ((m_dim==2 && isStructureManifold(n[0]/2, n[1]/2, newMaterialAssignment, mirroredStructure)) || (m_dim==3 && isStructureManifold(n[0]/2, n[1]/2, n[2]/2, newMaterialAssignment, N[0], N[1], N[2], mirroredStructure)) )
            {
              newParticuleGenerated = true;
            }
            else
            {
              newMaterialAssignment[indVoxel] = oldMat;
            }
          }
        }
      }
      if (newParticuleGenerated)
      {
        newMaterialAssignments.push_back(newMaterialAssignment);
        newVoxelsToProcess.push_back(voxelsToProcess[indParticule]);
      }
    }
    int istruct=0, nstruct=(int)newMaterialAssignments.size();
    for (istruct=0; istruct<nstruct; istruct++)
    {
      newMaterialAssignments[istruct] = getFullMaterialAssignment(n, newMaterialAssignments[istruct]);
    }
    
    // Compute Scores
    std::vector<cfgScalar> newParameters, newTensors;
    //std::vector<std::vector<std::vector<float> > > newx[3];
    //computeParameters(iStepperType, newMaterialAssignments, n, baseMaterialStructures[prevLevel], N, newParameters, newx); 
    computeParametersAndTensorValues(n, newMaterialAssignments, newParameters, newTensors);

    int nbPrevStructures = (int)allMaterialAssignments.size();
    allPhysicalParameters.insert(allPhysicalParameters.end(), newParameters.begin(), newParameters.end());
    allMaterialAssignments.insert(allMaterialAssignments.end(), newMaterialAssignments.begin(), newMaterialAssignments.end());
    allTensors.insert(allTensors.end(), newTensors.begin(), newTensors.end());

    //writeFiles(iLevel, materialAssignments, baseMaterialStructures[prevLevel], physicalParameters, x, "SMC_before_resampling"+std::to_string(icycle));
    writeFiles(iLevel, newMaterialAssignments, baseMaterialStructures[prevLevel], newParameters, newTensors, "SMC_new_samples_"+std::to_string(icycle));

    std::vector<cfgScalar> scores = scoringFunction.computeScores(allPhysicalParameters);
    scores = getSubVector(scores, genIncrementalSequence(nbPrevStructures, (int)scores.size()-1));

    // Resampling
    std::vector<int> newparticules;
    resampler.resample(scores, nTargetParticules, newparticules);
    nparticule = (int)newparticules.size();
    voxelsToProcess = getSubVector(newVoxelsToProcess, newparticules);
    particulesToProcess.clear();
    for (int iparticule=0; iparticule<nparticule; iparticule++)
    {
      if (voxelsToProcess[iparticule].size()>0)
      {
        particulesToProcess.insert(iparticule);
      }
    }
    physicalParameters = getSubVector(newParameters, dimparam, newparticules);
    materialAssignments = getSubVector(newMaterialAssignments, newparticules);
    tensors = getSubVector(newTensors, m_dim==2?6:21, newparticules);

    //writeFiles(iLevel, materialAssignments, baseMaterialStructures[prevLevel], physicalParameters, x, "SMC_after_resampling_"+std::to_string(icycle));
    writeFiles(iLevel, materialAssignments, baseMaterialStructures[prevLevel], physicalParameters, tensors, "SMC_after_resampling_"+std::to_string(icycle));

    icycle++;
  }
  oNewMaterialAssignments = allMaterialAssignments;
  oNewPhysicalParameters = allPhysicalParameters;
  oNewTensors = allTensors;

  std::cout<< "oNewPhysicalParameters" << oNewPhysicalParameters.size() << std::endl;
}


void SamplesGeneratorImpl::getNewCombinations(const std::vector<int> &iBoundaryPoints, const std::vector<cfgScalar> &iPoints, 
                                              int N[3], const std::vector<std::vector<int> > &iBaseMaterialStructures, int iNbCombinations, 
                                              std::vector<std::vector<int> > &oNewMaterialAssignments, std::vector<int> &oNewBaseMaterialStructures)
{
  std::vector<cfgScalar> boundaryPointCoords = getSubVector(iPoints, 3, iBoundaryPoints);
  cfgScalar bary[3];
  getBarycenter<cfgScalar,3>(boundaryPointCoords, bary);
  Vector3S c(bary[0],bary[1],bary[2]);

  cfgScalar SqMaxDist = 0;
  int ipoint=0, npoint=(int)iBoundaryPoints.size();
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    int indPoint = iBoundaryPoints[ipoint];
    Vector3S p = getVector3S(indPoint, iPoints);
    cfgScalar SqDist = (p-c).squaredNorm();
    if (SqDist>SqMaxDist)
    {
      SqMaxDist = SqDist;
    }
  }

  cfgScalar SqMaxDist2 = 0;
  npoint = iPoints.size()/3;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    Vector3S p = getVector3S(ipoint, iPoints);
    int closestBoundaryPoint = getClosestPoint(p, boundaryPointCoords);
    Vector3S q = getVector3S(iBoundaryPoints[closestBoundaryPoint], iPoints);
    float SqDist = (p-q).squaredNorm();
    if (SqDist>SqMaxDist2)
    {
      SqMaxDist2 = SqDist;
    }
  }
  SqMaxDist = SqMaxDist2;

  int maxBoundaryMismatch = 3*N[0]*N[0]*4*2;

  oNewBaseMaterialStructures.clear();
  oNewMaterialAssignments.clear();
  int nMat = iPoints.size()/3;
  int nVox = pow(2, m_dim);
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
        Vector3S p = getVector3S(matIndex, iPoints);
        int closestBoundaryPoint = getClosestPoint(p, boundaryPointCoords);
        Vector3S q = getVector3S(iBoundaryPoints[closestBoundaryPoint], iPoints);
        float SqDist = (p-q).squaredNorm();
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
      //std::cout << "# assignments = " << oNewMaterialAssignments.size() << std::endl;
    }
  }
  oNewBaseMaterialStructures = toStdVector(BaseMaterials);
}

void SamplesGeneratorImpl::getNewCombinationsV2(const std::vector<int> &iBoundaryPoints, const std::vector<cfgScalar> &iPoints, int N[3], int iNbCombinations, 
                                                std::vector<std::vector<int> > &oNewMaterialAssignments, std::vector<int> &oNewBaseMaterialStructures)
{
  std::vector<cfgScalar> boundaryPointCoords = getSubVector(iPoints, 3, iBoundaryPoints);

  oNewBaseMaterialStructures.clear();
  oNewMaterialAssignments.clear();

  int nMat = iPoints.size()/3;
  int nVox = pow(2, m_dim);
  std::vector<int> materialAssignement(nVox, -1);
  float pmin = 0.1;
  float pmax = 1;
  std::set<int> BaseMaterials;

  int ipoint=0, npoint=(int)iBoundaryPoints.size();
  cfgScalar SqMaxDist2 = 0;
  npoint = iPoints.size()/3;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    Vector3S p = getVector3S(ipoint, iPoints);
    int closestBoundaryPoint = getClosestPoint(p, boundaryPointCoords);
    Vector3S q = getVector3S(iBoundaryPoints[closestBoundaryPoint], iPoints);
    float SqDist = (p-q).squaredNorm();
    if (SqDist>SqMaxDist2)
    {
      SqMaxDist2 = SqDist;
    }
  }
  cfgScalar SqMaxDist = SqMaxDist2;

  std::vector<int> voxelIndices;
  int ivox;
  for (ivox=0; ivox<nVox; ivox++)
  {
    voxelIndices.push_back(ivox);
  }

  int ivertex=0, nvertex=(int)iBoundaryPoints.size();
  for (ivertex=0; ivertex<nvertex; ivertex++)
  {
    int pointIndex = iBoundaryPoints[ivertex];
    materialAssignement.clear();
    materialAssignement.resize(nVox, pointIndex);
    BaseMaterials.insert(pointIndex);
    oNewMaterialAssignments.push_back(materialAssignement);

    std::vector<int> materialAssignementInit = materialAssignement;

    int ivox=0;
    for (ivox=0; ivox<nVox; ivox++)
    {
      materialAssignement = materialAssignementInit;

      int indVox = ivox;

      int itrial=0, ntrial=100;
      for (itrial=0; itrial<ntrial; itrial++)
      {
        bool assignmentOk = true;
        bool matAssigned = false;  
        while (!matAssigned && assignmentOk)
        {
          int matIndex = rand() % nMat;
          Vector3S p = getVector3S(matIndex, iPoints);
          int closestBoundaryPoint = getClosestPoint(p, boundaryPointCoords);
          Vector3S q = getVector3S(iBoundaryPoints[closestBoundaryPoint], iPoints);
          float SqDist = (p-q).squaredNorm();
          float a = sqrt(SqDist/SqMaxDist);
          float proba = a*pmin + (1-a)*pmax;
          if (SqDist < 1.e-6)
          {
            proba = 1;
          }
          else
          {
            proba = 0.1;
          }
          float r = (float)rand()/float(RAND_MAX);
          if (proba > r)
          {
            materialAssignement[indVox] = matIndex;
            BaseMaterials.insert(matIndex);
            matAssigned = true;
          }
        }
        if (assignmentOk)
        {
          oNewMaterialAssignments.push_back(materialAssignement);
        }
      }
    }
  }
  oNewBaseMaterialStructures = toStdVector(BaseMaterials);
}

void SamplesGeneratorImpl::getNewCombinationsV3(const std::vector<int> &iBoundaryIndexArray, const std::vector<cfgScalar> &iPoints, int N[3], int iNbCombinations, 
                                              std::vector<std::vector<int> > &oNewMaterialAssignments, std::vector<int> &oNewBaseMaterialStructures)
{
  oNewBaseMaterialStructures.clear();
  oNewMaterialAssignments.clear();

  std::vector<cfgScalar> sampledBoundary;
  sampleMesh(iPoints, iBoundaryIndexArray, 1000, sampledBoundary);

  std::vector<cfgScalar> distToBoundary;
  computeDistancesToPointCloud(iPoints, sampledBoundary, distToBoundary);
  cfgScalar maxDist = *std::max_element(distToBoundary.begin(), distToBoundary.end()); 
  distToBoundary = cfgUtil::mult<cfgScalar>(distToBoundary, 1./maxDist);

  std::vector<cfgScalar> probabilities;
  int ipoint=0, npoint=(int)iPoints.size()/3;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    cfgScalar proba = 1-distToBoundary[ipoint];
    probabilities.push_back(proba);
  }

  int nMat = iPoints.size()/3;
  int nVox = pow(2, m_dim);
  std::vector<int> materialAssignement(nVox, -1);
  float pmin = 0.;
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
        float proba = probabilities[matIndex];

        float r = (float)rand()/float(RAND_MAX);
        if (proba > r)
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
      //std::cout << "# assignments = " << oNewMaterialAssignments.size() << std::endl;
    }
  }
  oNewBaseMaterialStructures = toStdVector(BaseMaterials);
}

int SamplesGeneratorImpl::getRandomMatWithLinearProba(std::multimap<float, int> &iDist2MatIndex, float eps, float r)
{
  std::multimap<float, int> & dist2Mat = iDist2MatIndex;

  float dist = (float)rand()/float(RAND_MAX);
  cfgScalar proba = 1-dist;
  if (proba > r)
  {
    std::vector<int> closestMat;
    std::multimap<float,int>::iterator it,itlow,itup;
    itlow = dist2Mat.lower_bound(dist); 
    itup = dist2Mat.upper_bound(dist); 

    if (itlow!=dist2Mat.begin())
    {
      itlow--;
    }
    float distlow = itlow->first;
    float distup = itup->first;
    bool same_value = sqrt(distup-distlow) < eps*eps;
    if (dist-distlow < distup-dist || same_value)
    {
      closestMat.push_back(itlow->second);
      dist = itlow->first;

      itlow--;
      for (it=itlow; ; it--)
      {
        if (it->first > dist-eps)
        {
          closestMat.push_back(it->second);
        }
        else
        {
          break;
        }
        if (it==dist2Mat.begin())
        {
          break;
        }
      }
    }
    if (dist-distlow >= distup-dist || same_value)
    {
      dist = itup->first;
      if (closestMat.size()==0 || itup->second != closestMat.front())
      {
        closestMat.push_back(itup->second);
      }
      itup++;
      for (it=itup; it!=dist2Mat.end(); it++)
      {
        if (it->first < dist+eps)
        {
          closestMat.push_back(it->second);
        }
        else
        {
          break;
        }
      }
    }
    int matIndex = closestMat[rand() % closestMat.size()];
    return matIndex;
  }
  else
  {
    return -1;
  }
}

void SamplesGeneratorImpl::getNewCombinationsV4(const std::vector<int> &iBoundaryIndexArray, const std::vector<cfgScalar> &iPoints, int N[3], int iNbCombinations, 
                                              std::vector<std::vector<int> > &oNewMaterialAssignments, const std::vector<std::vector<int> >  &iNewBaseMaterialStructures)
{
  oNewMaterialAssignments.clear();

  std::vector<cfgScalar> sampledBoundary;
  sampleMesh(iPoints, iBoundaryIndexArray, 1000, sampledBoundary);

  std::vector<cfgScalar> boundaryPoints = getSubVector(iPoints, 3, iBoundaryIndexArray);
  sampledBoundary.insert(sampledBoundary.end(), boundaryPoints.begin(), boundaryPoints.end());

  std::vector<cfgScalar> distToBoundary;
  computeDistancesToPointCloud(iPoints, sampledBoundary, distToBoundary);
  cfgScalar maxDist = *std::max_element(distToBoundary.begin(), distToBoundary.end()); 
  distToBoundary = cfgUtil::mult<cfgScalar>(distToBoundary, 1./maxDist);

  std::multimap<float, int> dist2Mat;
  int ipoint=0, npoint=(int)iPoints.size()/3;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    dist2Mat.insert(std::pair<float,int>(distToBoundary[ipoint],ipoint));
  }
  float eps = 1.e-4;

  int nMat = iPoints.size()/3;
  int nVox = pow(2, m_dim);
  std::vector<int> materialAssignement(nVox, -1);
  float pmin = 0.;
  float pmax = 1;
  std::set<int> BaseMaterials;
  int n[2] = {2*N[0], 2*N[1]};

  std::vector<int> voxelIndices;
  int ivox;
  for (ivox=0; ivox<nVox; ivox++)
  {
    voxelIndices.push_back(ivox);
  }
  bool mirroredStructure = m_orthotropicOnly||m_cubicOnly;
  while (oNewMaterialAssignments.size()<iNbCombinations)
  {
    materialAssignement.clear();
    materialAssignement.resize(nVox, -1);

    random_shuffle(voxelIndices.begin(), voxelIndices.end());

    bool assignmentOk = true;
    int ivox=0;
    for (ivox=0; ivox<nVox && assignmentOk; ivox++)
    {
      int indVox = voxelIndices[ivox];

      bool matAssigned = false;  
      while (!matAssigned && assignmentOk)
      {
        float r = (float)rand()/float(RAND_MAX);
        int matIndex = getRandomMatWithLinearProba(dist2Mat, eps, r);
        if (matIndex >= 0)
        {
          materialAssignement[indVox] = matIndex;
          matAssigned = true;
          BaseMaterials.insert(matIndex);
        }
      } 
    }
    if (assignmentOk)
    {
      std::vector<int> cellMaterials;
      getMaterialAssignment(2, 2, materialAssignement, N[0], N[1], iNewBaseMaterialStructures, m_blockRep, m_blockRep, m_nbSubdivisions, cellMaterials);

      if (isStructureManifold(n[0], n[1], cellMaterials, N[0], N[1], 4, mirroredStructure))
      {
        oNewMaterialAssignments.push_back(materialAssignement);
      }
      //oNewMaterialAssignments.push_back(materialAssignement);
      //std::cout << "# assignments = " << oNewMaterialAssignments.size() << std::endl;
    }
  } 
}

void SamplesGeneratorImpl::getNewCombinationsV5(const std::vector<int> &iBoundaryIndexArray, const std::vector<cfgScalar> &iPoints, int N[3], int iNbCombinations, 
                                              std::vector<std::vector<int> > &oNewMaterialAssignments, const std::vector<std::vector<int> >  &iNewBaseMaterialStructures)
{
  oNewMaterialAssignments.clear();

  std::vector<cfgScalar> sampledBoundary;
  sampleMesh(iPoints, iBoundaryIndexArray, 1000, sampledBoundary);

  std::vector<cfgScalar> boundaryPoints = getSubVector(iPoints, 3, iBoundaryIndexArray);
  sampledBoundary.insert(sampledBoundary.end(), boundaryPoints.begin(), boundaryPoints.end());

  std::vector<cfgScalar> distToBoundary;
  computeDistancesToPointCloud(iPoints, sampledBoundary, distToBoundary);
  cfgScalar maxDist = *std::max_element(distToBoundary.begin(), distToBoundary.end()); 
  distToBoundary = cfgUtil::mult<cfgScalar>(distToBoundary, 1./maxDist);

  std::multimap<float, int> dist2Mat;
  int ipoint=0, npoint=(int)iPoints.size()/3;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    dist2Mat.insert(std::pair<float,int>(distToBoundary[ipoint],ipoint));
  }
  float eps = 1.e-4;

  std::vector<int> allBoundaryPoints;
  std::multimap<float, int>::iterator it, it_end=dist2Mat.end();
  for (it=dist2Mat.begin(); it!=it_end; it++)
  {
    if (it->first < eps)
      allBoundaryPoints.push_back(it->second);
  }
  int nboundary = (int)allBoundaryPoints.size();

  int nMat = iPoints.size()/3;
  int nVox = pow(2, m_dim);
  std::vector<int> materialAssignement(nVox, -1);
  float pmin = 0.;
  float pmax = 1;
  std::set<int> BaseMaterials;
  int n[2] = {2*N[0], 2*N[1]};

  std::vector<int> voxelIndices;
  int ivox;
  for (ivox=0; ivox<nVox; ivox++)
  {
    voxelIndices.push_back(ivox);
  }
  bool mirroredStructure = m_orthotropicOnly||m_cubicOnly;
  int indBoundary = 0;
  while (oNewMaterialAssignments.size()<iNbCombinations)
  {
    int matIndex = allBoundaryPoints[indBoundary++];
    indBoundary = indBoundary%nboundary;

    materialAssignement.clear();
    materialAssignement.resize(nVox, matIndex);

    oNewMaterialAssignments.push_back(materialAssignement);
    BaseMaterials.insert(matIndex);

    random_shuffle(voxelIndices.begin(), voxelIndices.end());
    int ivox=0;
    for (ivox=0; ivox<nVox; ivox++)
    {
      int indVox = voxelIndices[ivox];

      bool matAssigned = false;  
      while (!matAssigned)
      {
        float dist = (float)rand()/float(RAND_MAX);
        cfgScalar proba = 1-dist;

        float r = (float)rand()/float(RAND_MAX);
        if (proba > r)
        {
          std::vector<int> closestMat;
          std::multimap<float,int>::iterator it,itlow,itup;
          itlow = dist2Mat.lower_bound(dist); 
          itup = dist2Mat.upper_bound(dist); 

          if (itlow!=dist2Mat.begin())
          {
            itlow--;
          }
          float distlow = itlow->first;
          float distup = itup->first;
          bool same_value = sqrt(distup-distlow) < eps*eps;
          if (dist-distlow < distup-dist || same_value)
          {
            closestMat.push_back(itlow->second);
            dist = itlow->first;

            itlow--;
            for (it=itlow; ; it--)
            {
              if (it->first > dist-eps)
              {
                closestMat.push_back(it->second);
              }
              else
              {
                break;
              }
              if (it==dist2Mat.begin())
              {
                break;
              }
            }
          }
          if (dist-distlow >= distup-dist || same_value)
          {
            dist = itup->first;
            if (closestMat.size()==0 || itup->second != closestMat.front())
            {
              closestMat.push_back(itup->second);
            }
            itup++;
            for (it=itup; it!=dist2Mat.end(); it++)
            {
              if (it->first < dist+eps)
              {
                closestMat.push_back(it->second);
              }
              else
              {
                break;
              }
            }
          }
          int matIndex = closestMat[rand() % closestMat.size()];
  
          materialAssignement[indVox] = matIndex;
          std::vector<int> cellMaterials;
          getMaterialAssignment(2, 2, materialAssignement, N[0], N[1], iNewBaseMaterialStructures, m_blockRep, m_blockRep, m_nbSubdivisions, cellMaterials);

          if (isStructureManifold(n[0], n[1], cellMaterials, N[0], N[1], 4, mirroredStructure))
          {
            oNewMaterialAssignments.push_back(materialAssignement);
            BaseMaterials.insert(matIndex);
            matAssigned = true;
          } 
        }
      } 
    }
  } 
}

int SamplesGeneratorImpl::getClosestPoint(Vector3S &iP, const std::vector<cfgScalar> &iPoints)
{
  int ClosestPointIndex = -1;
  float SqDistMin = FLT_MAX;
  int ipoint=0, npoint=(int)iPoints.size()/3;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    Vector3S Q = getVector3S(ipoint, iPoints);
    float SqDist = (Q-iP).squaredNorm();
    if (SqDist < SqDistMin)
    {
      SqDistMin = SqDist;
      ClosestPointIndex = ipoint;
    }
  }
  return ClosestPointIndex;
}

void SamplesGeneratorImpl::computeParameters(std::string & iStepperType, const std::vector<std::vector<int> > &iMaterialAssignments, int n[3], const std::vector<std::vector<int> > &iBaseMaterialAssignments, int N[3],
                                             std::vector<cfgScalar> &oParameters, std::vector<std::vector<std::vector<float> > > oX[3])
{
  std::vector<std::vector<float> > stresses[3];

  computeDeformation(iStepperType, iMaterialAssignments, n, iBaseMaterialAssignments, N, stresses, oX);

  int ncell[3] = {m_blockRep*n[0]*m_nbSubdivisions, m_blockRep*n[1]*m_nbSubdivisions, m_blockRep*n[2]*m_nbSubdivisions};
  
  if (m_dim==2)
  {
    std::vector<std::vector<float> > strains[2][2];
    cfgMaterialUtilities::computeStrain(ncell, oX, strains);
    cfgMaterialUtilities::computeMaterialParameters(iMaterialAssignments, iBaseMaterialAssignments, stresses, strains, oParameters);
  }
  else
  {
    std::vector<std::vector<float> > strains[3][3];
    cfgMaterialUtilities::computeStrain3D(ncell, oX, strains);
    cfgMaterialUtilities::computeMaterialParameters(iMaterialAssignments, iBaseMaterialAssignments, stresses, strains, oParameters);
  }
}

void SamplesGeneratorImpl::computeParametersAndTensorValues(int n[3], const std::vector<std::vector<int> > &iMatAssignments, std::vector<cfgScalar> &oParameters, std::vector<cfgScalar> &oTensorValues)
{
  oParameters.clear();
  oTensorValues.clear();

  int imat, nmat=(int)iMatAssignments.size();
  for (imat=0; imat<nmat; imat++)
  {
    std::cout << "************ imat = " << imat << std::endl;
    NumericalCoarsening numCoarsening;
    NumericalCoarsening::StructureType type = (m_cubicOnly? NumericalCoarsening::Cubic: NumericalCoarsening::General);
    if (m_dim==2)
    {
      numCoarsening.computeCoarsenedElasticityTensorAndParameters(iMatAssignments[imat], n, m_blockRep, m_nbSubdivisions, m_mat2D, type, oTensorValues, oParameters);
    }
    else
    {
      numCoarsening.computeCoarsenedElasticityTensorAndParameters(iMatAssignments[imat], n, m_blockRep, m_nbSubdivisions, m_mat, type, oTensorValues, oParameters);
    }
  }
}

void SamplesGeneratorImpl::computeParametersAndTensorValues(int n[3], const std::vector<std::vector<double> > &iMatDistributions, std::vector<cfgScalar> &oParameters, std::vector<cfgScalar> &oTensorValues)
{
  oParameters.clear();
  oTensorValues.clear();

  int ncell = (int)iMatDistributions[0].size();
  std::vector<int> matAssignment = genIncrementalSequence(0, ncell-1);

  int nmat=(int)iMatDistributions.size();
  for (int imat=0; imat<nmat; imat++)
  {
    if (imat % 100==0 )
      std::cout << "************ imat = " << imat << std::endl;
    NumericalCoarsening numCoarsening;
    NumericalCoarsening::StructureType type = NumericalCoarsening::General;
    if (m_cubicOnly)
    {
      type = NumericalCoarsening::Cubic;
    }
    else if (m_orthotropicOnly)
    {
      type = NumericalCoarsening::Orthotropic;
    }
 
    std::vector<MaterialQuad2D> mat2D(ncell);
    std::vector<StrainLin2D> ene(ncell);

    for (unsigned int ii = 0; ii < mat2D.size(); ii++)
    {
      double val = iMatDistributions[imat][ii];
      float coeff = (cfgScalar) (val / (3 - 2 * val));
      ene[ii].param[0] = coeff * m_mat2D[1].e[0]->param[0];
      ene[ii].param[1] = coeff * m_mat2D[1].e[0]->param[1];

      for (unsigned int jj = 0; jj < mat2D[ii].e.size(); jj++)
      {
        mat2D[ii].e[jj] = &ene[ii];
      }
    }
    std::vector<cfgScalar> densities = convertVec<double, cfgScalar>(iMatDistributions[imat]);
    if (m_dim==2)
    {
      numCoarsening.computeCoarsenedElasticityTensorAndParameters(matAssignment, n, m_blockRep, m_nbSubdivisions, mat2D, type, densities, oTensorValues, oParameters);
    }
    else
    {
      numCoarsening.computeCoarsenedElasticityTensorAndParameters(matAssignment, n, m_blockRep, m_nbSubdivisions, mat2D, type, oTensorValues, oParameters);
    }

    if (imat % 10==0 )
    {
      writeFiles(n[0], convertVec<double, cfgScalar>(iMatDistributions), oParameters, oTensorValues, "continuous_"+std::to_string(imat));
    }
  }
}


void SamplesGeneratorImpl::writeFiles(int iLevel, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<std::vector<int> > &iBaseMaterialStructures, const std::vector<cfgScalar> &iParameters,
                                      const std::vector<std::vector<std::vector<float> > > iX[2], const std::string iPostFix)
{
  std::string fileRootName = m_OutputDirectory + "level" + std::to_string(iLevel) + "_";
  if (iPostFix!="")
  {
    fileRootName += iPostFix + "_";
  }
  std::string fileExtension = ".bin";
  cfgUtil::writeBinary<float>(fileRootName + "params" + fileExtension, iParameters);
  cfgUtil::writeBinary<int>(fileRootName + "baseMat" + fileExtension, iBaseMaterialStructures);
  cfgUtil::writeBinary<int>(fileRootName + "matAssignments" + fileExtension, iMaterialAssignments);
}

void SamplesGeneratorImpl::writeFiles(int iLevel, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<std::vector<int> > &iBaseMaterialStructures, const std::vector<cfgScalar> &iParameters,
                                      const std::vector<cfgScalar> &iTensorsValues, const std::string iPostFix)
{
  std::string fileRootName = m_OutputDirectory + "level" + std::to_string(iLevel) + "_";
  if (iPostFix!="")
  {
    fileRootName += iPostFix + "_";
  }
  std::string fileExtension = ".bin";
  cfgUtil::writeBinary<float>(fileRootName + "params" + fileExtension, iParameters);
  cfgUtil::writeBinary<int>(fileRootName + "baseMat" + fileExtension, iBaseMaterialStructures);
  cfgUtil::writeBinary<int>(fileRootName + "matAssignments" + fileExtension, iMaterialAssignments);
  cfgUtil::writeBinary<float>(fileRootName + "elasticityTensors" + fileExtension, iTensorsValues);
}

void SamplesGeneratorImpl::writeFiles(int iLevel, const std::vector<std::vector<float> > &iMaterialDistributions, const std::vector<cfgScalar> &iParameters, const std::vector<cfgScalar> &iTensorsValues, const std::string iPostFix)
{
  std::string fileRootName = m_OutputDirectory + "level" + std::to_string(iLevel) + "_";
  if (iPostFix!="")
  {
    fileRootName += iPostFix + "_";
  }
  std::string fileExtension = ".bin";
  cfgUtil::writeBinary<float>(fileRootName + "params" + fileExtension, iParameters);
  cfgUtil::writeBinary<float>(fileRootName + "matDistributions" + fileExtension, iMaterialDistributions);
  cfgUtil::writeBinary<float>(fileRootName + "elasticityTensors" + fileExtension, iTensorsValues);
}


bool SamplesGeneratorImpl::readFiles(int iLevel, std::vector<std::vector<int> > &oMaterialAssignments, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<cfgScalar> &oParameters)
{
  oMaterialAssignments.clear();
  oBaseMaterialStructures.clear();
  oParameters.clear();

  int version = 1;
  bool ResOk = true;
  if (version==0)
  {
    int readBinary = true;
    int blockSize = m_blockRep;
    std::string materialFile = m_OutputDirectory + "Materials_" + std::to_string(iLevel) + ".txt";
    std::string baseFileName = "StressDeformation_" + std::to_string(iLevel) + "_" + std::to_string(blockSize);
    std::string extension = (readBinary? ".bin": ".txt");
    std::string stressDeformationFileNames[2] = {m_OutputDirectory+ baseFileName + "_x" + extension, m_OutputDirectory + baseFileName + "_y" + extension};
    ResOk = cfgMaterialUtilities::computeMaterialParametersFromDeformations(materialFile, stressDeformationFileNames, readBinary, oParameters, oBaseMaterialStructures, oMaterialAssignments);
  }
  else
  {
    std::string fileRootName = m_OutputDirectory + "level" + std::to_string(iLevel) + "_";
    std::string fileExtension = ".bin";

    ResOk = cfgUtil::readBinary<float>(fileRootName + "params" + fileExtension, oParameters);
    if (ResOk)
      ResOk = cfgUtil::readBinary<int>(fileRootName + "baseMat" + fileExtension, oBaseMaterialStructures);
    if (ResOk)
      ResOk = cfgUtil::readBinary<int>(fileRootName + "matAssignments" + fileExtension, oMaterialAssignments);
  }
  return ResOk;
}

bool SamplesGeneratorImpl::readFiles(int iLevel, std::vector<std::vector<int> > &oMaterialAssignments, std::vector<std::vector<int> > &oBaseMaterialStructures, std::vector<cfgScalar> &oParameters, std::vector<cfgScalar> &oTensors,
                                     const std::string iPostFix)
{
  oMaterialAssignments.clear();
  oBaseMaterialStructures.clear();
  oParameters.clear();
  oTensors.clear();

  bool ResOk = true;

  std::string fileRootName = m_OutputDirectory + "level" + std::to_string(iLevel) + "_";
  if (iPostFix!="")
  {
    fileRootName += iPostFix + "_";
  }
  std::string fileExtension = ".bin";

  ResOk = cfgUtil::readBinary<float>(fileRootName + "params" + fileExtension, oParameters);
  if (ResOk)
    ResOk = cfgUtil::readBinary<int>(fileRootName + "baseMat" + fileExtension, oBaseMaterialStructures);
  if (ResOk)
    ResOk = cfgUtil::readBinary<int>(fileRootName + "matAssignments" + fileExtension, oMaterialAssignments);
  if (ResOk)
    ResOk = cfgUtil::readBinary<float>(fileRootName + "elasticityTensors" + fileExtension, oTensors);

  return ResOk;
}

void SamplesGeneratorImpl::concatenateFiles(int iLevel, int indexMin, int indexMax, const std::string iPostFix, const std::string iNewPostFix)
{
  bool ResOk = true;

  std::vector<std::vector<int> > allMatAssignments;
  std::vector<std::vector<int> > allBaseMaterialStructures;
  std::vector<cfgScalar> allParameters;
  std::vector<cfgScalar> allTensors;

  for (int i=indexMin; i<=indexMax; i++)
  {
    std::vector<std::vector<int> > matAssignments;
    std::vector<std::vector<int> > baseMaterialStructures;
    std::vector<cfgScalar> parameters;
    std::vector<cfgScalar> tensors;
    ResOk = readFiles(iLevel, matAssignments, baseMaterialStructures, parameters, tensors, iPostFix + std::to_string(i));
    if (ResOk)
    {
      allMatAssignments.insert(allMatAssignments.end(), matAssignments.begin(), matAssignments.end());
      allBaseMaterialStructures.insert(allBaseMaterialStructures.end(), baseMaterialStructures.begin(), baseMaterialStructures.end());
      allParameters.insert(allParameters.end(), parameters.begin(), parameters.end());
      allTensors.insert(allTensors.end(), tensors.begin(), tensors.end());
    }
    else
    {
      std::cout << "failed to read files for index = " << i << std::endl;
    }
  }
  writeFiles(iLevel, allMatAssignments, allBaseMaterialStructures, allParameters, allTensors, iNewPostFix);
}

int SamplesGeneratorImpl::run()
{
  if (0)
  {
    int level = 16;
    int indexMin = 0;
    int indexMax = 19;
    std::string postFix = "SMC_after_resampling";
    std::string newPostFix = "";
    concatenateFiles(level, indexMin, indexMax, postFix, newPostFix);
    exit(0);
  }
  if (0)
  {
    std::string inputFileName =  m_OutputDirectory + "x_level2//StressStrain_2_3";
    std::string outputFileName =  m_OutputDirectory + "x_level2//StressStrain_2_3_x.txt";
    concatenateData(inputFileName, 65536, outputFileName);
    exit(0);
  }
  if (0)
  {
    std::string inputFileName = m_OutputDirectory + "x_level1//StressStrain_1_3_x.txt";
    Vector3S forceAxis;
    std::vector<std::vector<int> > materials;
    std::vector<std::vector<float> > stresses, strains;
    bool ResOk = readData(inputFileName,  forceAxis, materials, stresses, strains);
    exit(0);
  }

  srand(0);

  float E1, nu1, lambda1, mu1;
  mu1 = 100;
  lambda1 = 1000;
  fromLamesParametersToYoungModulusPoissonRatio(lambda1, mu1, E1, nu1);
 
  float E2, nu2, lambda2, mu2;
  E2 = E1/1000;
  nu2 = 0.1;
  fromYoungModulusPoissonRatioToLamesParameters(E2, nu2, lambda2, mu2);

  //std::vector<StrainEneNeo> ene(2);
  std::vector<StrainLin2D> ene2D(2);
  /*ene2D[0].param[0] = 1;
  ene2D[0].param[1] = 10;
  ene2D[1].param[0] = 10;
  ene2D[1].param[1] = 100;*/ 
  /*ene2D[0].param[0] = 0.1;
  ene2D[0].param[1] = 1;
  ene2D[1].param[0] = 100;
  ene2D[1].param[1] = 1000;*/ 
  ene2D[0].param[0] = mu2;
  ene2D[0].param[1] = lambda2;
  ene2D[1].param[0] = mu1;
  ene2D[1].param[1] = lambda1;

  //std::vector<StrainEneNeo> ene(2);
  std::vector<StrainLin> ene(2);
  ene[0].param[0] = mu2;
  ene[0].param[1] = lambda2;
  ene[1].param[0] = mu1;
  ene[1].param[1] = lambda1;

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

  // mirror structures
  if (0)
  {
    int level = 1;
    for (int ilevel=1; ilevel<=16; ilevel*=2)
    {
      level = ilevel;
      std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
      std::vector<float> physicalParameters;
      bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters);

      int n[2] = {level, level};
      int nmat = (int)materialAssignments.size();
      std::vector<std::vector<int> > newMaterialAssignments(nmat);
      for (int imat=0; imat<nmat; imat++)
      {
        mirrorStructure(n[0], n[1], materialAssignments[imat], newMaterialAssignments[imat]);
      }
      std::string fileRootName = m_OutputDirectory + "level" + std::to_string(2*level) + "_";
      std::string fileExtension = ".bin";
      //fileRootName += "new_";
      cfgUtil::writeBinary<int>(fileRootName + "matAssignments" + fileExtension, newMaterialAssignments);
    }
    return 0;
  }


  // compute material properties
  if (0)
  {
    int level = 32;
    std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
    std::vector<float> physicalParameters;
    bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters);

    int N[2] = {1,1};
    int n[2] = {level, level};
    std::vector<cfgScalar> newParameters;
    std::vector<std::vector<std::vector<float> > > newx[2];
    computeParameters(stepperType, materialAssignments, n, baseMaterialStructures, N, newParameters, newx); 
    
    std::string fileRootName = m_OutputDirectory + "level" + std::to_string(level) + "_";
    std::string fileExtension = ".bin";
    //fileRootName += "new_";
    cfgUtil::writeBinary<float>(fileRootName + "params" + fileExtension, newParameters);
    return 0;
  }

  // compute homogenised tensors
  if (0)
  {
    int level = 4;
    std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
    std::vector<float> physicalParameters;
    bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters);

    int n[3] = {level, level, level};
    std::vector<float> tensors;
    std::vector<float> newphysicalParameters;
    computeParametersAndTensorValues(n, materialAssignments, newphysicalParameters, tensors);
    writeFiles(level, materialAssignments, baseMaterialStructures, newphysicalParameters, tensors);
    return 0;
  }

  // compute homogenised tensors from continuous material distribution
  if (0)
  {
    int level = 16;

    std::vector<std::vector<double> > materialDistributions;
    bool ResOk = cfgUtil::readBinary<double>(m_OutputDirectory + "16Ortho.bin", materialDistributions);
    if (ResOk)
    {
      int n[3] = {level, level, level};
      std::vector<float> tensors;
      std::vector<float> physicalParameters;
      computeParametersAndTensorValues(n, materialDistributions, physicalParameters, tensors);
      writeFiles(level, convertVec<double, cfgScalar>(materialDistributions), physicalParameters, tensors, "continuous");
    }
    return 0;
  }

  //computeVariations(stepperType, 10);

  std::vector<std::vector<int> > baseMaterialStructures(2);
  baseMaterialStructures[0].push_back(0);
  baseMaterialStructures[1].push_back(1);
  int N[3] = {1,1,1};

  if (0)
  {
    int maxlevel = 2;
    if (m_orthotropicOnly)
      maxlevel = 4;
    if (m_cubicOnly)
      maxlevel = 8; //maxlevel = (m_dim==2? 8: 4);

    for (int ilevel=1; ilevel<=maxlevel; ilevel*=2)
    {
      int n[3] = {ilevel, ilevel, ilevel};
      std::vector<std::vector<int> > newMaterialAssignment;
      computeMaterialParametersIncremental(stepperType, ilevel, newMaterialAssignment);
      //computeMaterialParameters(stepperType, newMaterialAssignment, n, baseMaterialStructures, N, ilevel, m_blockRep, true, true);

      std::vector<float> tensors;
      std::vector<float> physicalParameters;
      //std::vector<std::vector<std::vector<float> > > x[3];
      //computeParameters(stepperType, newMaterialAssignment, n, baseMaterialStructures, N, physicalParameters, x); 
      computeParametersAndTensorValues(n, newMaterialAssignment, physicalParameters, tensors);
      writeFiles(ilevel, newMaterialAssignment, baseMaterialStructures, physicalParameters, tensors);
    }
  }

  // Continuous optimization
  if (0)
  {
    int ilevel = 16;
    int n[3] = {ilevel, ilevel, ilevel};

    std::vector<std::vector<int> > materialAssignments, baseMaterialStructures, clusters;
    std::vector<float> physicalParameters, tensors;

    bool ResOk = readFiles(ilevel, materialAssignments, baseMaterialStructures, physicalParameters, tensors);
    std::vector<std::vector<int> > newMaterialAssignments;
    std::vector<cfgScalar> newParameters, newTensors;
    runContinuousOptimization(ilevel, m_mat2D, materialAssignments, physicalParameters, tensors, newMaterialAssignments, newParameters, newTensors);
  }

  if (1)
  {
    bool growAllStructures = false;

    for (int ilevel=16; ilevel<=16; ilevel *=2)
      //for (int ilevel=4; ilevel<=16; ilevel*=2)
        //for (int ilevel=16; ilevel<=32; ilevel*=2)
          //for (int ilevel=4; ilevel<=128; ilevel*=2)
    {
      int n[3] = {ilevel, ilevel, ilevel};

      std::vector<std::vector<std::vector<int> > > materialAssignments, baseMaterialStructures, clusters;
      std::vector<std::vector<float> > physicalParameters, tensors;
      int nlevel=ilevel;
      materialAssignments.resize(nlevel+1);
      baseMaterialStructures.resize(nlevel+1);
      physicalParameters.resize(nlevel+1);
      tensors.resize(nlevel+1);

      // Read previous files
      // ===================
      for (int iprevlevel=1; iprevlevel<nlevel; iprevlevel*=2)
      {
        //bool ResOk = readFiles(iprevlevel, materialAssignments[iprevlevel], baseMaterialStructures[iprevlevel], physicalParameters[iprevlevel]);
        bool ResOk = readFiles(iprevlevel, materialAssignments[iprevlevel], baseMaterialStructures[iprevlevel], physicalParameters[iprevlevel], tensors[iprevlevel]);
      }

      std::vector<std::vector<int> > allNewMaterialAssignments;
      std::vector<cfgScalar> allParameters, allTensors;
      std::vector<std::vector<std::vector<float> > > x[2];
      int prevLevel = ilevel/2;

      bool useSMC = true;
      if (useSMC)
      {
        computeVariationsSMC(ilevel, stepperType, materialAssignments, baseMaterialStructures, physicalParameters, tensors, allNewMaterialAssignments, allParameters, allTensors); 
        baseMaterialStructures[ilevel] = baseMaterialStructures[prevLevel];
      }
      else
      {
        // Grow structures
        // ===============
        if (growAllStructures)
        {
          std::cout << "Grow structures" << std::endl;
          int n_prev[3] = {prevLevel, prevLevel, prevLevel};  
          growStructureDoubleSize(n_prev,  materialAssignments[prevLevel], materialAssignments[ilevel]);
          std::cout << "nb init comb = " << materialAssignments[prevLevel].size() <<" nb comb = " << materialAssignments[ilevel].size() << std::endl;
          baseMaterialStructures[ilevel] = baseMaterialStructures[prevLevel];
          computeParameters(stepperType, materialAssignments[ilevel], n, baseMaterialStructures[ilevel], N, physicalParameters[ilevel], x); 
          if (1)
          {
            writeFiles(ilevel, materialAssignments[ilevel], baseMaterialStructures[prevLevel], physicalParameters[ilevel], x, "scaled");
          }
        }
        else
        {
          baseMaterialStructures[ilevel] = baseMaterialStructures[prevLevel];
        }

        // Compute variations
        // ==================
        std::cout << "Compute variations" << std::endl;
        int maxLevel=ilevel;
        for (int isublevel=1; isublevel<maxLevel; isublevel*=2)
        {
          std::vector<std::vector<int> > newMatAssignments;
          if (growAllStructures)
          {
            computeVariationsV2(ilevel, isublevel, materialAssignments, baseMaterialStructures, physicalParameters, newMatAssignments);
          }
          else
          {
            computeVariationsV3(ilevel, isublevel, materialAssignments, baseMaterialStructures, physicalParameters, newMatAssignments);
          }

          std::cout << "Nb new assignments for level " << ilevel << " and sublevel " << isublevel << " = " << newMatAssignments.size() << std::endl;

          std::vector<cfgScalar> newParameters;
          std::vector<std::vector<std::vector<float> > > newx[2];
          computeParameters(stepperType, newMatAssignments, n, baseMaterialStructures[ilevel], N, newParameters, newx); 

          writeFiles(ilevel, newMatAssignments, baseMaterialStructures[ilevel], newParameters, x, "sublevel"+ std::to_string(isublevel));
          allNewMaterialAssignments.insert(allNewMaterialAssignments.end(), newMatAssignments.begin(), newMatAssignments.end());
          allParameters.insert(allParameters.end(), newParameters.begin(), newParameters.end());
          x[0].insert(x[0].end(), newx[0].begin(), newx[0].end());
          x[1].insert(x[1].end(), newx[1].begin(), newx[1].end());
        } 
      }
      physicalParameters[ilevel] = allParameters;

      //computeParameters(stepperType, allNewMaterialAssignments, n, baseMaterialStructures[ilevel], N, physicalParameters[ilevel], x); 
      //writeFiles(ilevel, allNewMaterialAssignments, baseMaterialStructures[ilevel], physicalParameters[ilevel], x, "non_clustered");
      writeFiles(ilevel, allNewMaterialAssignments, baseMaterialStructures[ilevel], physicalParameters[ilevel], allTensors, "non_clustered");

      //writeFiles(ilevel, newMaterialAssignment, baseMaterialStructures, physicalParameters, tensors);

      // Cluster
      std::cout << "Reduce" << std::endl;
      bool reduce = false;
      if (reduce && !useSMC)
      {
        int nbClusters = 2000;
        std::vector<cfgScalar> lengths(3, 1);
        int nbIter = 10;

        std::vector<float> rescaledParameters = physicalParameters[ilevel];
        rescaleData(rescaledParameters, 3, lengths);

        std::vector<std::vector<int> >  clusters;
        std::vector<int> pointIndices; 
        getKMeans(nbIter, nbClusters, rescaledParameters, 3, clusters, &pointIndices);

        for (int iaxis=0; iaxis<2; iaxis++)
        {
          x[iaxis] = getSubVector(x[iaxis], pointIndices);
        } 
        physicalParameters[ilevel] = getSubVector(physicalParameters[ilevel], 3, pointIndices);
        allNewMaterialAssignments = getSubVector(allNewMaterialAssignments, pointIndices);
      }
      //writeStressDeformationFile(allNewMaterialAssignments, n, stresses, x);
      writeFiles(ilevel, allNewMaterialAssignments, baseMaterialStructures[ilevel], physicalParameters[ilevel], x);
    } 
  }
 
  bool writeSingleFile = true, readSingleFile = true;
  bool writeFullDeformation = true, readFullDeformation = true;
  if (0)
  {
    //computeMaterialParametersLevel1(stepperType, writeSingleFile, writeFullDeformation);
    //computeMaterialParameters(stepperType, 2, (m_dim==2?16:256), readSingleFile, writeSingleFile, readFullDeformation, writeFullDeformation);
    //computeMaterialParameters(stepperType, 3, 20);
    
    //computeMaterialParameters(stepperType, 3, 65536, readSingleFile, writeSingleFile, readFullDeformation, writeFullDeformation);

    //computeMaterialParameters(material, stepperType, 3, 586);
    //computeMaterialParameters(material, stepperType, 3, 429);

    //computeMaterialParametersIncremental(stepperType, 4);
    //computeMaterialParametersIncremental(stepperType, 8);
    /*computeMaterialParametersIncremental(stepperType, 1);
    computeMaterialParametersIncremental(stepperType, 2);
    computeMaterialParametersIncremental(stepperType, 3);*/ 
    int ilevel;
    for (ilevel=3; ilevel<7; ilevel++)
    {
      //int indLevel = ilevel;
      int indLevel = pow(2, ilevel);
      //computeMaterialParametersIncremental(stepperType, indLevel);
    }

    return 0 ;
  }
  return 0;
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

        float forceMagnitude = 1.f;
        int NumberOfSample = 10;
        int axis = 0;
        Vector3S forceAxis;
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




