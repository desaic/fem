
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

#include "SamplerUtilities.h"
using namespace SamplerUtilities;

#include "FamilyGenerator.h"

#include "ScoringFunction.h"
#include "Resampler.h"
#include "NumericalCoarsening.h"
#include "NumericalCoarseningUsingMultigrid.h"
#include "MaterialOptimizer.h"
#include "DistanceField.h"

#include "ThermalAnalysis.h"
#include "StrengthAnalysis.h"

#include "ConfigFile.hpp"

enum Type
{
  Cubic2D,
  Orthotropic2D,
  Cubic3D
};

enum Stage
{
  ExhaustiveGamutComputationStage,
  DiscreteAndContinuousOptimizationStage,
  DiscreteOptimizationStage,
  ContinuousOptimizationStage,
  StrengthAnalysisStage,
  VariantsComputationStage,
  FilesConcatenationStage,
  IsotropicToOrthotropicConversionStage,
  SubsamplingStage,
  MirrorStructureStage,
  From2Dto3DConversionStage,
  UpscaleMicrostructureStage,
  FixDisconnectedMicrostructureStage,
  FixNonManifoldMicrostructureStage,
  MaterialPropertiesComputationStage,
  TensorHomogenizationStage,
  FamilyGenerationStage,
  ThermalAnalysisStage,
};

SamplesGeneratorImpl::SamplesGeneratorImpl()
{
  m_UseLinearMaterial = true;

  m_blockRep = 1;// (m_dim==2? 8: 4);
  m_nbSubdivisions = 1;
  m_continuousMatDist = false;

  m_orthotropicOnly = false;
  m_cubicOnly = false;
  m_filterOutDisconnetedStructures = true;
  m_useLogScale = true;

  m_useVariableDensity = true;

  m_nmat = 2;
  m_dim = 2;
  setOutputDirectory("..//..//Output_tmp//");
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
    newBaseMaterialStructures.resize(m_nmat);
    for (int imat=0; imat<m_nmat; imat++)
    {
      std::vector<int> mat(1, imat);
      newBaseMaterialStructures[imat] = mat;
    }
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
      int nmat = m_nmat;
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
          if (isStructureManifold(nprev[0], nprev[1], nprev[2], matAssignmentMirroredAlongDiag, N[0], N[1], N[2], true, m_nmat))
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
      int nmat = m_nmat;
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
          if (isStructureManifold(nprev[0], nprev[1], nprev[2], matAssignment, N[0], N[1], N[2], true, m_nmat))
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
            if (isStructureManifold(n[0], n[1], n[2], matAssignment, N[0], N[1], N[2], false, m_nmat))
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
  
  std::set<std::vector<int> > setNewMaterials;
  int indVoidMat = std::min_element(m_baseMaterialDensities.begin(), m_baseMaterialDensities.end()) - m_baseMaterialDensities.begin();

  for (int imat=0; imat<(int)newMaterialAssignments.size(); imat++)
  {
    if (imat % 10000 == 1)
    {
      std::cout << "imat = " << imat << std::endl;
    }
    std::vector<int> & newMaterialAssignment = newMaterialAssignments[imat];
    if (m_filterOutDisconnetedStructures)
    {
      std::vector<int> binaryMaterialAssignment = toBinaryDistribution(newMaterialAssignment, indVoidMat);
      bool resOk = true;
      bool modified=false;
      if (m_dim==2)
      {
        resOk = filterOutNonConnectedComponents(n[0], n[1], binaryMaterialAssignment, modified);
      }
      else
      {
        resOk = filterOutNonConnectedComponents(n[0], n[1], n[2], binaryMaterialAssignment, modified);
      }
      if (resOk && !modified)
      {
        setNewMaterials.insert(newMaterialAssignment);
      }
    }
    else
    {
      setNewMaterials.insert(newMaterialAssignment);
    }
  }
  newMaterialAssignments = toStdVector(setNewMaterials);
  std::cout << "** Nb manifold comb = " << newMaterialAssignments.size() << std::endl;

  oNewMaterialAssignments = newMaterialAssignments;

  return ResOk;
}

std::vector<int> SamplesGeneratorImpl::toBinaryDistribution(const std::vector<int> &iMatAssignment, int iVoidMatIndex)
{
  std::vector<int> binaryMatAssignment;

  std::vector<int> materials(m_nmat, 1);
  materials[iVoidMatIndex] = 0;
  int ivoxel=0, nvoxel=(int)iMatAssignment.size();
  for (ivoxel=0; ivoxel<nvoxel; ivoxel++)
  {
    int indMat = iMatAssignment[ivoxel];
    int newIndMat = materials[indMat];
    binaryMatAssignment.push_back(newIndMat);
  }
  return binaryMatAssignment;
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
      std::vector<int> matAssignmentQuarter;
      getQuarter(n[0], n[1], iMatAssignments, matAssignmentQuarter);
      getTriangularStructure(n[0]/2, n[1]/2, matAssignmentQuarter, matAssignments);
    }
    else
    {
      std::vector<int> matAssignmentQuarter;
      getQuarter(n[0], n[1], n[2], iMatAssignments, matAssignmentQuarter);
      getTetrahedralStructure(n[0]/2, n[1]/2, n[2]/2, matAssignmentQuarter, matAssignments);
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

int SamplesGeneratorImpl::runContinuousOptimization(int iLevel, int iStartCycle, const std::string &iSuffix)
{
  int level = iLevel;
  int startCycle = iStartCycle;

  bool fixNonManifoldStructure = true;

  std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
  std::vector<float> physicalParameters, tensors;
  bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters, tensors, iSuffix);

  std::vector<std::vector<int> > newMaterialAssignments;
  std::vector<cfgScalar> newParameters, newTensors;
  runContinuousOptimization(level, startCycle, fixNonManifoldStructure, materialAssignments, physicalParameters, tensors, newMaterialAssignments, newParameters, newTensors);
  writeFiles(level, newMaterialAssignments, baseMaterialStructures, newParameters, newTensors, "ContinuousOptimizationResult");

  /*
  std::vector<float>  targetParameters;
  bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters, tensors, "ContinuousOptimizationInit_0");
  ResOk = readFiles(level, materialAssignments, baseMaterialStructures, targetParameters, tensors, "ContinuousOptimizationTarget_0");

  std::vector<int> newparticules = genIncrementalSequence(0, 0);
  std::vector<std::vector<int> > newMaterialAssignments;
  std::vector<cfgScalar> newParameters, newTensors;
  runContinuousOptimization(level, startCycle, fixNonManifoldStructure, materialAssignments, physicalParameters, tensors, targetParameters, newparticules, newMaterialAssignments, newParameters, newTensors);
  */ 

  return 0;
}

void SamplesGeneratorImpl::runContinuousOptimization(int iLevel, int iCycle, bool iFixNonManifoldStructure,
                                                     const std::vector<std::vector<int> > &iMaterialAssignments, 
                                                     const std::vector<float> &iParameters,
                                                     const std::vector<float> &iTensors,
                                                     std::vector<std::vector<int> > &oNewMaterialAssignments,
                                                     std::vector<cfgScalar> &oNewPhysicalParameters,
                                                     std::vector<cfgScalar> &oNewTensors)
{
  int paramdim = getNbParameters();
  DistanceField distanceField(paramdim);
  std::vector<cfgScalar> derivatives;
  std::cout << "Computing distances..." << std::endl;
  std::vector<cfgScalar> distances = distanceField.computeDistances(iParameters, &derivatives);
  cfgScalar coef = 0.1;
  std::vector<cfgScalar> newPoints = cfgUtil::add<cfgScalar>(iParameters, cfgUtil::mult<cfgScalar>(derivatives, coef));
  

  cfgScalar minRadius = 0;
  int nTargetParticules = 100;
  minRadius = 0.1;
  if (!m_cubicOnly)
  {
    nTargetParticules = 1000;
    //minRadius = 0.05;
  }

  std::vector<int> newparticules;
  Resampler resampler;
  std::cout << "Resampling boundary..." << std::endl;
  resampler.resampleBoundary(minRadius, paramdim, iParameters, distances, nTargetParticules, newparticules);
  std::cout << "nb points = " << newparticules.size() << std::endl;

  /*int ind=0;
  std::cout << "Current point: ";
  for (int icoord=0; icoord<paramdim; icoord++)
  {
    std::cout <<  iParameters[icoord] << " ";
  }
  std::cout << std::endl;
  std::cout << "Target point: ";
  for (int icoord=0; icoord<paramdim; icoord++)
  {
    std::cout <<  newPoints[icoord] << " ";
  }
  std::cout << std::endl;*/ 

  std::vector<std::vector<int> > baseMaterialStructures(2);
  baseMaterialStructures[0].push_back(0);
  baseMaterialStructures[1].push_back(1);

  std::vector<cfgScalar> initParameters, initTensors;
  initParameters = getSubVector(iParameters, paramdim, newparticules);
  initTensors = getSubVector(iTensors, m_dim==2?6:21, newparticules);
  std::vector<std::vector<int> > &initMaterialAssignments = getSubVector(iMaterialAssignments, newparticules);
  writeFiles(iLevel, initMaterialAssignments, baseMaterialStructures, initParameters, initTensors, "ContinuousOptimizationInit_"+std::to_string(iCycle));

  std::vector<cfgScalar> targetParameters;
  targetParameters = getSubVector(newPoints, paramdim, newparticules);
  writeFiles(iLevel, initMaterialAssignments, baseMaterialStructures, targetParameters, initTensors, "ContinuousOptimizationTarget_"+std::to_string(iCycle));

  runContinuousOptimization(iLevel, iCycle, iFixNonManifoldStructure, iMaterialAssignments, iParameters, iTensors, newPoints, newparticules, oNewMaterialAssignments, oNewPhysicalParameters, oNewTensors);
}

void SamplesGeneratorImpl::runContinuousOptimization(int iLevel, int iCycle, bool iFixNonManifoldStructure,
                                                     const std::vector<std::vector<int> > &iMaterialAssignments, 
                                                     const std::vector<float> &iParameters,
                                                     const std::vector<float> &iTensors,
                                                     const std::vector<float> &iOffsetedParameters,
                                                     const std::vector<int> &iNewParticules,
                                                     std::vector<std::vector<int> > &oNewMaterialAssignments,
                                                     std::vector<cfgScalar> &oNewPhysicalParameters,
                                                     std::vector<cfgScalar> &oNewTensors)
{
  oNewMaterialAssignments.clear();
  oNewPhysicalParameters.clear();
  oNewTensors.clear();

  std::vector<std::vector<int> > allNewMaterialAssignments;
  std::vector<cfgScalar> allNewParameters;
  std::vector<cfgScalar> allNewTensors;

  const std::vector<int> &newparticules = iNewParticules;
  const std::vector<cfgScalar> &newPoints = iOffsetedParameters;

  int paramdim = getNbParameters();
  std::vector<std::vector<int> > baseMaterialStructures(2);
  baseMaterialStructures[0].push_back(0);
  baseMaterialStructures[1].push_back(1);

  MaterialOptimizer optimizer;
  if (m_cubicOnly)
  {
    optimizer.setStructureType(MaterialOptimizer::Cubic);
  }
  else if (m_orthotropicOnly)
  {
    optimizer.setStructureType(MaterialOptimizer::Orthotropic);
  }
  std::cout << "Running continuous optimization..." << std::endl;
  bool writeIntermediateResult = true;
  int N[3] = {iLevel, iLevel, (m_dim==2?1: iLevel)};
  int istruct=0, nstruct=(int)newparticules.size();
  //nstruct = 1;
  for (istruct=0; istruct<nstruct; istruct++)
  {
    std::cout << "mat #" << istruct << std::endl;
    int indStruct = newparticules[istruct];
    //int indStruct = 3;
    std::vector<int> matAssignment = iMaterialAssignments[indStruct];
    std::vector<float> targetParameters, currentParameters;
    for (int icoord=0; icoord<paramdim; icoord++)
    {
      targetParameters.push_back(newPoints[paramdim*indStruct+icoord]);
      currentParameters.push_back(iParameters[paramdim*indStruct+icoord]);
    }
    std::vector<std::vector<int> > newMaterials;
     std::vector<cfgScalar> newParameters, newTensors;
    bool resOk = true;
    if (m_dim==2)
    {
      resOk = optimizer.run2D(N, m_mat2D, matAssignment, currentParameters, targetParameters, newMaterials); 
    }
    else
    {
      resOk = optimizer.run3D(N, m_mat, matAssignment, currentParameters, targetParameters, newMaterials); 
    }
    newMaterials = getNewElements(allNewMaterialAssignments, newMaterials);
    if (resOk && newMaterials.size()>0)
    {
      if (iFixNonManifoldStructure)
      {
        int imat, nmat=(int)newMaterials.size();
        for (imat=0; imat<nmat; imat++)
        {
          fixNonManifoldStructure(N, newMaterials[imat]);
        }
      }
      allNewMaterialAssignments.insert(allNewMaterialAssignments.end(), newMaterials.begin(), newMaterials.end());
      if (writeIntermediateResult)
      {
        computeParametersAndTensorValues(N, newMaterials, newParameters, newTensors);
        allNewParameters.insert(allNewParameters.end(), newParameters.begin(), newParameters.end());
        allNewTensors.insert(allNewTensors.end(), newTensors.begin(), newTensors.end());
        writeFiles(iLevel, allNewMaterialAssignments, baseMaterialStructures, allNewParameters, allNewTensors, "ContinousOptimizationResult_tmp");
      }
    }
  }
  if (!writeIntermediateResult)
  {
    computeParametersAndTensorValues(N, allNewMaterialAssignments, allNewParameters, allNewTensors);
  }
  oNewMaterialAssignments = iMaterialAssignments;
  oNewPhysicalParameters = iParameters;
  oNewTensors = iTensors;

  oNewMaterialAssignments.insert(oNewMaterialAssignments.end(), allNewMaterialAssignments.begin(), allNewMaterialAssignments.end());
  oNewPhysicalParameters.insert(oNewPhysicalParameters.end(), allNewParameters.begin(), allNewParameters.end());
  oNewTensors.insert(oNewTensors.end(), allNewTensors.begin(), allNewTensors.end());
}

bool SamplesGeneratorImpl::fixNonManifoldStructure2D(int n[2], std::vector<int> &ioMatAssignments)
{
  bool fixed = false;

  if (0)
  {
    std::cout << "Init structure: " << std::endl;
    dumpStructure(n[0], n[1], ioMatAssignments);
  }
  bool mirroredStructure = m_orthotropicOnly||m_cubicOnly;

  std::vector<int> matAssignments;
  if (mirroredStructure)
  {
    getQuarter(n[0], n[1], ioMatAssignments, matAssignments);
  }
  int nx = n[0]/2;
  int ny = n[1]/2;

  std::vector<int> nonManifoldVertices;
  getNonManifoldVertices(nx, ny, matAssignments, mirroredStructure, nonManifoldVertices);
  if (nonManifoldVertices.size()>0)
  {
    std::set<int> setTriangleElements;
    if (m_cubicOnly)
    {
      std::vector<int> triangleElements;
      getTriangleElements(nx, ny, matAssignments, triangleElements);
      std::vector<int> nonManifoldVerticesTriangle;
      setTriangleElements = toStdSet(triangleElements);
      int ivertex, nvertex=(int)nonManifoldVertices.size();
      for (ivertex=0; ivertex<nvertex; ivertex++)
      {
        int indVertex = nonManifoldVertices[ivertex];
        if (setTriangleElements.count(indVertex)>0)
        {
          nonManifoldVerticesTriangle.push_back(indVertex);
        }
      }
      nonManifoldVertices = nonManifoldVerticesTriangle;
    }
    std::set<int> setNonManifoldVertices = toStdSet(nonManifoldVertices);
    while (setNonManifoldVertices.size()>0)
    {
      int indVertex = *setNonManifoldVertices.begin();
      setNonManifoldVertices.erase(indVertex);

      int x, y;
      getVectorIndexToGrid(indVertex, nx, ny, x, y);

      std::vector<int> cellMat0;
      int indMat11 = getGridToVectorIndex(x, y, nx, ny);
      int indMat12 = getGridToVectorIndex(x, (y+ny-1)%ny, nx, ny);
      int indMat21 = getGridToVectorIndex((x+nx-1)%nx, y, nx, ny);
      int indMat22 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, nx, ny);
      if (matAssignments[indMat11]==0)
        cellMat0.push_back(indMat11);
      if (matAssignments[indMat12]==0)
        cellMat0.push_back(indMat12);
      if (matAssignments[indMat21]==0)
        cellMat0.push_back(indMat21);
      if (matAssignments[indMat22]==0)
        cellMat0.push_back(indMat22);
      if (cellMat0.size()==2)
      {
        int index = rand()%2;
        int indVertexQuarter = cellMat0[index];
        if (m_cubicOnly && setTriangleElements.count(indVertexQuarter)==0)
        {
          indVertexQuarter = cellMat0[1-index];
        }
        matAssignments[indVertexQuarter] = 1;

        if (m_cubicOnly)
        {
          std::vector<int> matAssignmentTriangle;
          getTriangularStructure(nx, ny, matAssignments, matAssignmentTriangle);
          mirrorStructureAlongDiagonal(nx, ny, matAssignmentTriangle, matAssignments);
        }

        int x, y;
        getVectorIndexToGrid(indVertexQuarter, nx, ny, x, y);
        if (x>0 & y>0 && !isVertexManifold(x, y, nx, ny, matAssignments))
        {
          int newVertex = getGridToVectorIndex(x, y, nx, ny);
          setNonManifoldVertices.insert(newVertex);
        }
        if (x!=nx-1 && y>0 && !isVertexManifold((x+1)%nx, y, nx, ny, matAssignments))
        {
          int newVertex = getGridToVectorIndex((x+1)%nx, y, nx, ny);
          setNonManifoldVertices.insert(newVertex);
        }
        if (x>0 && y!=ny-1 && !isVertexManifold(x, (y+1)%ny, nx, ny, matAssignments))
        {
          int newVertex = getGridToVectorIndex(x, (y+1)%ny, nx, ny);
          setNonManifoldVertices.insert(newVertex);
        }
        if (x!=nx-1 && y!=ny-1 && !isVertexManifold((x+1)%nx, (y+1)%ny, nx, ny, matAssignments))
        {
          int newVertex = getGridToVectorIndex((x+1)%nx, (y+1)%ny, nx, ny);
          setNonManifoldVertices.insert(newVertex);
        }
      }
    }
    mirrorStructure(nx, ny, matAssignments, ioMatAssignments);
    if (0)
    {
      std::cout << "Result: " << std::endl;
      dumpStructure(n[0], n[1], ioMatAssignments);
    }
    fixed = true;
  }
  return fixed;
}

bool SamplesGeneratorImpl::fixNonManifoldStructure(int n[3], std::vector<int> &ioMatAssignments)
{
  if (m_dim==2)
    return fixNonManifoldStructure2D(n, ioMatAssignments);

  bool fixed = false;
  if (0)
  {
    std::cout << "Init structure: " << std::endl;
    dumpStructure(n[0], n[1], n[2], ioMatAssignments);
  }
  bool mirroredStructure = m_orthotropicOnly||m_cubicOnly;

  std::vector<int> matAssignments;
  if (mirroredStructure)
  {
    getQuarter(n[0], n[1], n[2], ioMatAssignments, matAssignments);
  }
  int nx = n[0]/2;
  int ny = n[1]/2;
  int nz = n[2]/2;

  std::vector<int> nonManifoldVertices;
  getNonManifoldVertices(nx, ny, nz, matAssignments, mirroredStructure, nonManifoldVertices);
  if (nonManifoldVertices.size()>0)
  {
    std::set<int> setTetElements;
    if (m_cubicOnly)
    {
      std::vector<int> tetElements;
      getTetrahedralElements(nx, ny, nz, matAssignments, tetElements);
      std::vector<int> nonManifoldVerticesTet;
      setTetElements = toStdSet(tetElements);
      int ivertex, nvertex=(int)nonManifoldVertices.size();
      for (ivertex=0; ivertex<nvertex; ivertex++)
      {
        int indVertex = nonManifoldVertices[ivertex];
        if (setTetElements.count(indVertex)>0)
        {
          nonManifoldVerticesTet.push_back(indVertex);
        }
      }
      nonManifoldVertices = nonManifoldVerticesTet;
    }
    std::set<int> setNonManifoldVertices = toStdSet(nonManifoldVertices);
    while (setNonManifoldVertices.size()>0)
    {
      int indVertex = *setNonManifoldVertices.begin();
      setNonManifoldVertices.erase(indVertex);

      std::vector<int> indMat;
      int x, y, z;
      getVectorIndexToGrid(indVertex, nx, ny, nz, x, y, z);
      if (!isVertexManifold(x, y, z, nx, ny, nz, matAssignments, mirroredStructure, &indMat))
      {
        std::vector<int> cellMat0;
        for (int imat=0; imat<indMat.size(); imat++)
        {
          if (matAssignments[indMat[imat]]==0)
          {
            cellMat0.push_back(indMat[imat]);
          }
        }
        int index = rand()%cellMat0.size();
        int indVertexQuarter = cellMat0[index];

        int x, y, z;
        getVectorIndexToGrid(indVertexQuarter, nx, ny, nz, x, y, z);

        if (m_cubicOnly && setTetElements.count(indVertexQuarter)==0)
        {
          if (z > y) std::swap(z,y);
          if (y > x) std::swap(y,x);
          if (z > y) std::swap(z,y);
          indVertexQuarter = getGridToVectorIndex(x, y, z, nx, ny, nz);
        }
        matAssignments[indVertexQuarter] = 1;

        if (m_cubicOnly)
        {
          std::vector<int> matAssignmentTet;
          getTetrahedralStructure(nx, ny, nz, matAssignments, matAssignmentTet);
          mirrorStructureAlongDiagonal(nx, ny, nz, matAssignmentTet, matAssignments);
        }

        int newVertex = getGridToVectorIndex(x, y, z, nx, ny, nz);
        setNonManifoldVertices.insert(newVertex);
        newVertex = getGridToVectorIndex(x, (y+1)%ny, z, nx, ny, nz);
        setNonManifoldVertices.insert(newVertex);
        newVertex = getGridToVectorIndex((x+1)%nx, y, z, nx, ny, nz);
        setNonManifoldVertices.insert(newVertex);
        newVertex = getGridToVectorIndex((x+1)%nx, (y+1)%ny, z, nx, ny, nz);
        setNonManifoldVertices.insert(newVertex);
        newVertex = getGridToVectorIndex(x, y, (z+1)%nz, nx, ny, nz);
        setNonManifoldVertices.insert(newVertex);
        newVertex = getGridToVectorIndex(x, (y+1)%ny, (z+1)%nz, nx, ny, nz);
        setNonManifoldVertices.insert(newVertex);
        newVertex = getGridToVectorIndex((x+1)%nx, y, (z+1)%nz, nx, ny, nz);
        setNonManifoldVertices.insert(newVertex);
        newVertex = getGridToVectorIndex((x+1)%nx, (y+1)%ny, (z+1)%nz, nx, ny, nz);
        setNonManifoldVertices.insert(newVertex);
      }
    }
    mirrorStructure(nx, ny, nz, matAssignments, ioMatAssignments);
    if (0)
    {
      std::cout << "Result: " << std::endl;
      dumpStructure(n[0], n[1], n[2], ioMatAssignments);
    }
    fixed = true;
  }
  return fixed;
}

void SamplesGeneratorImpl::computeVariationsSMC(int iLevel, 
                                               std::string & iStepperType,
                                               const std::vector<std::vector<int> > &iMaterialAssignments,  
                                               const std::vector<std::vector<int> >  &iBaseMaterialStructures, 
                                               const std::vector<float> &iPhysicalParameters,
                                               const std::vector<float> &iTensors,
                                               bool iGrowStructures,
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

  const std::vector<std::vector<int> > & baseMaterialStructures = iBaseMaterialStructures;

  std::vector<float> physicalParameters = iPhysicalParameters;
  std::vector<float> tensors = iTensors;

  int nTargetParticules = 1000;
  if (!m_cubicOnly)
  {
    nTargetParticules = 3000;
  }
  int dimparam = getNbParameters();
  int nInitParticules = (int)physicalParameters.size()/dimparam;

  // Sampling
  std::vector<int> particules;
  std::vector<cfgScalar> scores;

  bool useDistanceField = true;
  bool cacheDensities = true;
  ScoringFunction scoringFunction(dimparam);
  scoringFunction.setUseDistanceField(useDistanceField);
  scoringFunction.setCacheDensities(cacheDensities);
  scoringFunction.setUseLogScale(m_useLogScale);
  scoringFunction.setParametersToUse(m_parametersToUse);

  float minRadius = 0.;
  //float minRadius = 0.05;
  Resampler resampler;

  std::vector<std::vector<int> > materialAssignments, allMaterialAssignments;
  std::vector<float> allPhysicalParameters;
  std::vector<float> allTensors;
  std::set<std::vector<int> > materialSet;

  cfgScalar ratio = 0;

  bool filterOutDisconnectedStructures = m_filterOutDisconnetedStructures;
  bool growStructure = iGrowStructures;
  std::vector<cfgScalar> initParams;
  if (growStructure)
  {
    std::cout << "growing structures..." << std::endl;
      int istruct=0, nstruct = (int)iMaterialAssignments.size();
      for (istruct=0; istruct<nstruct; istruct++)
      {
        std::vector<int> materialAsssignmentPrevLevel = iMaterialAssignments[istruct];
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
        if (materialSet.count(initMaterialAssignment)==0)
        {
          materialSet.insert(initMaterialAssignment);
          allMaterialAssignments.push_back(initMaterialAssignment);
          for (int iparam=0; iparam<dimparam; iparam++)
          {
            initParams.push_back(physicalParameters[dimparam*istruct+iparam]);
          }
        }
      }
      bool recomputeProperties = true;
      if (recomputeProperties)
      {
        std::vector<cfgScalar> newParameters, newTensors;
        std::cout << "recomputing properties..." << std::endl;
        computeParametersAndTensorValues(n, allMaterialAssignments, newParameters, newTensors);
        physicalParameters = newParameters;
        writeFiles(iLevel, allMaterialAssignments, baseMaterialStructures, physicalParameters, tensors, "recomputed_properties");
      }
      else
      {
        physicalParameters = initParams;
      }
      std::cout << "resampling particules..." << std::endl;
      scoringFunction.setNewPointIndexStart(0);
      scores = scoringFunction.computeScores(physicalParameters);
      resampler.resample(minRadius, dimparam, physicalParameters, scores, nTargetParticules, particules, ratio);
  }
  else
  {
    materialSet = toStdSet(iMaterialAssignments);
    allMaterialAssignments = iMaterialAssignments;

    scoringFunction.setNewPointIndexStart(0);
    scores = scoringFunction.computeScores(physicalParameters);
    resampler.resample(minRadius, dimparam, physicalParameters, scores, nTargetParticules, particules, ratio);
  }
  materialAssignments = allMaterialAssignments;
  allPhysicalParameters = physicalParameters;
  allTensors = tensors;

  int nparticule = (int)particules.size();
  std::vector<std::vector<int> > voxelsToProcess(nparticule);

  int nvoxel = getNbFreeCells(n);
  std::vector<int> initVoxelIndices = genIncrementalSequence(0, nvoxel-1);

  //int nswap = nvoxel/8;
  int nswap = 1;
  std::cout << nswap << "/" << nvoxel << " to swap " << std::endl;

  for (int iparticule=0; iparticule<nparticule; iparticule++)
  {
    voxelsToProcess[iparticule] = initVoxelIndices;
  }
  physicalParameters = getSubVector(physicalParameters, dimparam, particules);
  tensors = getSubVector(tensors, m_dim==2?6:21, particules);
  materialAssignments = getSubVector(materialAssignments, 1, particules);

  writeFiles(iLevel, materialAssignments, baseMaterialStructures, physicalParameters, tensors, "SMC_init");

  bool mirroredStructure = m_orthotropicOnly||m_cubicOnly;

  std::cout << "processing particules..." << std::endl;
  std::set<int> particulesToProcess = toStdSet(genIncrementalSequence(0, nparticule-1));
  int icycle=0;
  while (particulesToProcess.size())
  {
    std::cout << "nb particules to process = " << particulesToProcess.size() << std::endl;
    std::vector<std::vector<int> > newVoxelsToProcess;
    std::vector<std::vector<int> > newMaterialAssignments;
    std::set<int>::iterator it, it_end=particulesToProcess.end();
    for (it=particulesToProcess.begin(); it!=it_end; it++)
    {
      int indParticule = *it;
 
      std::vector<int> newMaterialAssignment = getFreeCellsMaterialAssignment(n, materialAssignments[indParticule]);

      bool newParticuleGenerated = false;
      int indSwap = 0;
      while (voxelsToProcess[indParticule].size() && (!newParticuleGenerated || indSwap<nswap))
      {
        int randomIndex = rand() % voxelsToProcess[indParticule].size();
        int indVoxel = voxelsToProcess[indParticule][randomIndex];
        voxelsToProcess[indParticule].erase(voxelsToProcess[indParticule].begin()+randomIndex);

        int newMat = rand()%m_nmat;

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
            if ((m_dim==2 && isStructureManifold(n[0]/2, n[1]/2, mirroredMaterials, mirroredStructure)) || (m_dim==3 && isStructureManifold(n[0]/2, n[1]/2, n[2]/2, mirroredMaterials, N[0], N[1], N[2], mirroredStructure, m_nmat)) )
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
            if ((m_dim==2 && isStructureManifold(n[0]/2, n[1]/2, newMaterialAssignment, mirroredStructure)) || (m_dim==3 && isStructureManifold(n[0]/2, n[1]/2, n[2]/2, newMaterialAssignment, N[0], N[1], N[2], mirroredStructure, m_nmat)) )
            {
              newParticuleGenerated = true;
            }
            else
            {
              newMaterialAssignment[indVoxel] = oldMat;
            }
          }
        }
        indSwap++;
      }
      if (newParticuleGenerated)
      {
        newMaterialAssignments.push_back(newMaterialAssignment);
        newVoxelsToProcess.push_back(voxelsToProcess[indParticule]);
      }
    }
    int indVoidMat = std::min_element(m_baseMaterialDensities.begin(), m_baseMaterialDensities.end()) - m_baseMaterialDensities.begin();
    std::vector<std::vector<int> > newStructures;
    int istruct=0, nstruct=(int)newMaterialAssignments.size();
    for (istruct=0; istruct<nstruct; istruct++)
    {
      std::vector<int> newMaterialAssignment = getFullMaterialAssignment(n, newMaterialAssignments[istruct]);
      if (materialSet.count(newMaterialAssignment)==0)
      {
        bool resOk = true;
        if (filterOutDisconnectedStructures)
        {
          bool modified = false;
          if (m_dim==2)
          {
            std::vector<int> binaryMaterialAssignment = toBinaryDistribution(newMaterialAssignment, indVoidMat);
            resOk = filterOutNonConnectedComponents(n[0], n[1], binaryMaterialAssignment, modified);
            //resOk = filterOutNonConnectedComponents(n[0], n[1], newMaterialAssignment, modified);
          }
          else
          {
            std::vector<int> binaryMaterialAssignment = toBinaryDistribution(newMaterialAssignment, indVoidMat);
            resOk = filterOutNonConnectedComponents(n[0], n[1], n[2], binaryMaterialAssignment, modified);
            //resOk = filterOutNonConnectedComponents(n[0], n[1], n[2], newMaterialAssignment, modified);
          }
          if (modified && materialSet.count(newMaterialAssignment)>0 )
          {
            resOk = false;
          }
        }
        if (resOk)
        {
          newStructures.push_back(newMaterialAssignment);
          materialSet.insert(newMaterialAssignment);
        }
      } 
    }
    newMaterialAssignments = newStructures; 
    
    // Compute Scores
    std::vector<cfgScalar> newParameters, newTensors;
    //std::vector<std::vector<std::vector<float> > > newx[3];
    //computeParameters(iStepperType, newMaterialAssignments, n, baseMaterialStructures[prevLevel], N, newParameters, newx); 
    computeParametersAndTensorValues(n, newMaterialAssignments, newParameters, newTensors);

    int nbPrevStructures = (int)allMaterialAssignments.size();
    allPhysicalParameters.insert(allPhysicalParameters.end(), newParameters.begin(), newParameters.end());
    allMaterialAssignments.insert(allMaterialAssignments.end(), newMaterialAssignments.begin(), newMaterialAssignments.end());
    allTensors.insert(allTensors.end(), newTensors.begin(), newTensors.end());

    writeFiles(iLevel, newMaterialAssignments, baseMaterialStructures, newParameters, newTensors, "SMC_new_samples_"+std::to_string(icycle));

    std::set<std::vector<int> > setMaterials = toStdSet(newMaterialAssignments);
    std::cout << "Nb different particules = " << setMaterials.size() << std::endl;
    if (setMaterials.size() < 10)
    {
      break;
    }

    scoringFunction.setNewPointIndexStart(nbPrevStructures);
    std::vector<cfgScalar> scores = scoringFunction.computeScores(allPhysicalParameters);
    scores = getSubVector(scores, genIncrementalSequence(nbPrevStructures, (int)scores.size()-1));

    // Resampling
    std::vector<int> newparticules;
    resampler.resample(scores, nTargetParticules, newparticules);
    //resampler.resample(minRadius, dimparam, newParameters, scores, nTargetParticules, newparticules);
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

    //writeFiles(iLevel, materialAssignments, baseMaterialStructures, physicalParameters, tensors, "SMC_after_resampling_"+std::to_string(icycle));
    writeFiles(iLevel, allMaterialAssignments, baseMaterialStructures, allPhysicalParameters, allTensors, "SMC_Result_tmp");

    icycle++;
  }
  oNewMaterialAssignments = allMaterialAssignments;
  oNewPhysicalParameters = allPhysicalParameters;
  oNewTensors = allTensors;

  //std::cout<< "oNewPhysicalParameters =" << oNewPhysicalParameters.size() << std::endl;
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

  if (m_dim==2)
  {
    NumericalCoarsening::StructureType type = (m_cubicOnly? NumericalCoarsening::Cubic: NumericalCoarsening::General);
    NumericalCoarsening numCoarsening;
    numCoarsening.init(n, m_blockRep, m_nbSubdivisions, m_mat2D);
    int imat, nmat=(int)iMatAssignments.size();
    for (imat=0; imat<nmat; imat++)
    {
      if (imat % 100 == 0)
        std::cout << "************ imat = " << imat << std::endl;
      numCoarsening.computeCoarsenedElasticityTensorAndParameters2D(iMatAssignments[imat], n, m_baseMaterialDensities, type, oTensorValues, oParameters);
    }
  }
  else
  {
    bool useMultigridSolver = true;
    if (useMultigridSolver)
    {
      NumericalCoarseningUsingMultigrid::StructureType type = (m_cubicOnly? NumericalCoarseningUsingMultigrid::Cubic: NumericalCoarseningUsingMultigrid::General);
      NumericalCoarseningUsingMultigrid numCoarsening[8];
      for (int i=0; i<8; i++)
      {
        numCoarsening[i].setUseVariableDensity(m_useVariableDensity);
        numCoarsening[i].init(n, m_blockRep, m_nbSubdivisions, m_mat, type);
      }
      int paramdim = getNbParameters();
      int imat, nmat=(int)iMatAssignments.size();
      oParameters.resize(nmat*paramdim);
      for (imat=0; imat<nmat; imat+=8)
      {
#pragma omp parallel for
        for (int i=0; i<8; i++)
        {
          int indMat = imat+i;
          if (indMat < nmat)
          {
            //std::cout << "************ imat = " << indMat << " max param = " << (oParameters.size()>0? *std::max_element(oParameters.begin(), oParameters.end()): 0) << " min param = "<< (oParameters.size()>0? *std::min_element(oParameters.begin(), oParameters.end()): 0) << std::endl;
            std::cout << "************ imat = " << indMat << std::endl;
            std::vector<cfgScalar> parameters;
            numCoarsening[i].computeCoarsenedElasticityTensorAndParameters3D(iMatAssignments[indMat], n, m_baseMaterialDensities, type, oTensorValues, parameters);
            for (int iparam=0; iparam<paramdim; iparam++)
            {
              std::cout << parameters[iparam] << " ";
              oParameters[indMat*paramdim + iparam] = parameters[iparam];
            } 
            std::cout << std::endl;
          }
        }
      }
    }
    else
    {
      NumericalCoarsening::StructureType type = (m_cubicOnly? NumericalCoarsening::Cubic: NumericalCoarsening::General);
      NumericalCoarsening numCoarsening;
      numCoarsening.init(n, m_blockRep, m_nbSubdivisions, m_mat);
      int imat, nmat=(int)iMatAssignments.size();
      for (imat=0; imat<nmat; imat++)
      {
        if (imat % 10 == 0)
          std::cout << "************ imat = " << imat << std::endl;
        numCoarsening.computeCoarsenedElasticityTensorAndParameters3D(iMatAssignments[imat], n, m_baseMaterialDensities, type, oTensorValues, oParameters);
      }
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

  bool readTextFiles = true;
  if (readTextFiles)
  {
    std::string fileRootName = m_OutputDirectory + "level" + std::to_string(iLevel) + "_";
    std::string fileExtension = ".txt";
    std::string fileName = fileRootName + "matAssignments" + fileExtension;
    std::ifstream stream(fileName);
    bool resOk = stream.is_open();
    if (resOk)
    {
      int nstruct = 13587;
      //int nstruct = 13689;
      int nvoxel = 32*32;
      int nvalues = nvoxel*nstruct;
      std::vector<int> mat(nvalues);
      deSerialize<int>(stream, nvalues, &mat[0]);

      int ind=0;
      oMaterialAssignments.reserve(nstruct);
      for (int istruct=0; istruct<nstruct; istruct++)
      {
        int sumMat = 0;
        std::vector<int> matAssignment;
        for (int ivoxel=0; ivoxel<nvoxel; ivoxel++)
        {
          int matIndex = mat[ind++];
          sumMat += matIndex;
          matAssignment.push_back(matIndex);
        }
        //if (sumMat>0)
        {
          oMaterialAssignments.push_back(matAssignment);
        }
      }
      std::cout << "nb valid structures = " << oMaterialAssignments.size() << std::endl;
    }
    return resOk;
  }

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

  ResOk = cfgUtil::readBinary<int>(fileRootName + "matAssignments" + fileExtension, oMaterialAssignments);
  if (ResOk)
    ResOk = cfgUtil::readBinary<int>(fileRootName + "baseMat" + fileExtension, oBaseMaterialStructures);
  if (ResOk)
    ResOk = cfgUtil::readBinary<float>(fileRootName + "params" + fileExtension, oParameters);
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

void SamplesGeneratorImpl::concatenateFiles(int iLevel, const std::string iPostFix1, const std::string iPostFix2, const std::string iNewPostFix)
{
  bool ResOk = true;

  std::vector<std::vector<int> > allMatAssignments;
  std::vector<std::vector<int> > allBaseMaterialStructures;
  std::vector<cfgScalar> allParameters;
  std::vector<cfgScalar> allTensors;

  ResOk = readFiles(iLevel, allMatAssignments, allBaseMaterialStructures, allParameters, allTensors, iPostFix1);
  if (ResOk)
  {
    std::vector<std::vector<int> > matAssignments;
    std::vector<std::vector<int> > baseMaterialStructures;
    std::vector<cfgScalar> parameters;
    std::vector<cfgScalar> tensors;
    ResOk = readFiles(iLevel, matAssignments, baseMaterialStructures, parameters, tensors, iPostFix2);
    if (ResOk)
    {
      allMatAssignments.insert(allMatAssignments.end(), matAssignments.begin(), matAssignments.end());
      allBaseMaterialStructures.insert(allBaseMaterialStructures.end(), baseMaterialStructures.begin(), baseMaterialStructures.end());
      allParameters.insert(allParameters.end(), parameters.begin(), parameters.end());
      allTensors.insert(allTensors.end(), tensors.begin(), tensors.end());
    }
    else
    {
      std::cout << "failed to read files 2"<< std::endl;
    }
    writeFiles(iLevel, allMatAssignments, allBaseMaterialStructures, allParameters, allTensors, iNewPostFix);
  }
  else
  {
    std::cout << "failed to read files 1"<< std::endl;
  }
}

bool SamplesGeneratorImpl::readDisneyFiles(int iDim, bool iCubic, bool iOrthotropic, std::vector<cfgScalar> &oPhysicalParameters)
{
  bool resOk = true;

  if (iDim==2 && iCubic)
  {
    resOk = resOk && readDisneyFiles("2D_cubic_o", iDim, iCubic, iOrthotropic, oPhysicalParameters);
    resOk = resOk && readDisneyFiles("2D_cubic_diamond", iDim, iCubic, iOrthotropic, oPhysicalParameters);
    resOk = resOk && readDisneyFiles("2D_cubic_cube", iDim, iCubic, iOrthotropic, oPhysicalParameters);
    resOk = resOk && readDisneyFiles("2D_cubic_rhombus", iDim, iCubic, iOrthotropic, oPhysicalParameters);
    resOk = resOk && readDisneyFiles("2D_cubic_rhomb-o", iDim, iCubic, iOrthotropic, oPhysicalParameters);
  }
  else if (iDim==2 && iOrthotropic)
  {
    resOk = resOk && readDisneyFiles("2D_ortho_o", iDim, iCubic, iOrthotropic, oPhysicalParameters);
    resOk = resOk && readDisneyFiles("2D_ortho_x", iDim, iCubic, iOrthotropic, oPhysicalParameters);
    resOk = resOk && readDisneyFiles("2D_ortho_hourglass", iDim, iCubic, iOrthotropic, oPhysicalParameters);
    resOk = resOk && readDisneyFiles("2D_ortho_diamond", iDim, iCubic, iOrthotropic, oPhysicalParameters);
  }
  else if (iDim==3 && iCubic)
  {
    resOk = resOk && readDisneyFiles("3D_cubic_x", iDim, iCubic, iOrthotropic, oPhysicalParameters);
    resOk = resOk && readDisneyFiles("3D_cubic_cross", iDim, iCubic, iOrthotropic, oPhysicalParameters);
    resOk = resOk && readDisneyFiles("3D_cubic_strut", iDim, iCubic, iOrthotropic, oPhysicalParameters);
  }
  return resOk;
}

bool SamplesGeneratorImpl::readDisneyFiles(const std::string &iFileName, int iDim, bool iCubic, bool iOrthotropic, std::vector<cfgScalar> &ioPhysicalParameters)
{
  std::string fileName = m_OutputDirectory + iFileName + ".txt";

  std::ifstream stream(fileName);
  if (!stream.is_open())
  {
    return false;
  }

  int nbStruct = 0;
  deSerialize<int>(stream, 1, &nbStruct);

  int dimparam = 0;
  if (iCubic)
    dimparam = 4;
  else if (iOrthotropic && iDim==2)
    dimparam = 5;

  int nvalues = nbStruct*dimparam;
  std::vector<float> values(nvalues);
  deSerialize<float>(stream, nvalues, &values[0]);

  if (iCubic)
  {
    //Relative Young's modulus,	Poisson's ratio,	Relative shear modulus,	Density
    assert(values.size()%dimparam==0);
    int ival, nval=(int)values.size();
    for (ival=0; ival<nval; ival+=dimparam)
    {
      cfgScalar E= values[ival];
      cfgScalar nu = values[ival+1];
      cfgScalar mu = values[ival+2];
      cfgScalar density = values[ival+3];

      ioPhysicalParameters.push_back(density);
      ioPhysicalParameters.push_back(E);
      ioPhysicalParameters.push_back(nu);
      ioPhysicalParameters.push_back(mu);
    }
  }
  else if (iOrthotropic && iDim==2)
  {
    // Relative Young's modulus X,	Relative Young's modulus Y,	Poisson's ratio YX,	Relative shear modulus,	Density
    assert(values.size()%dimparam==0);
    int ival, nval=(int)values.size();
    for (ival=0; ival<nval; ival+=dimparam)
    {
      cfgScalar E_x= values[ival];
      cfgScalar E_y= values[ival+1];
      if (E_x>1.e-6 && E_y>1.e-6)
      {
        cfgScalar nu_yx = values[ival+2];
        cfgScalar mu_xy = values[ival+3];
        cfgScalar density = values[ival+4];
        cfgScalar nu_xy = (E_y>1.e-6? nu_yx*E_x/E_y: 0);

        ioPhysicalParameters.push_back(density);
        ioPhysicalParameters.push_back(E_x);
        ioPhysicalParameters.push_back(E_y);
        ioPhysicalParameters.push_back(nu_xy);
        ioPhysicalParameters.push_back(mu_xy);
      }
    }
  }
  else
  {
    // not implemented
    assert(0);
  }
  return true;
}

bool SamplesGeneratorImpl::readPanettaFiles(int iDim, std::vector<cfgScalar> &oPhysicalParameters)
{
  bool resOk = true;

  if (iDim==2)
  {
    resOk = readPanettaFiles("moduli_2d", iDim, oPhysicalParameters);
  }
  else if (iDim==3)
  {
    resOk = readPanettaFiles("moduli_3d", iDim, oPhysicalParameters);
  }
  return resOk;
}

bool SamplesGeneratorImpl::readPanettaFiles(const std::string &iFileName, int iDim, std::vector<cfgScalar> &ioPhysicalParameters)
{
  std::string fileName = m_OutputDirectory + iFileName + ".txt";

  std::ifstream stream(fileName);
  if (!stream.is_open())
  {
    return false;
  }

  int nbStruct = 0;
  deSerialize<int>(stream, 1, &nbStruct);

  int dimparam = 3;

  int nvalues = nbStruct*dimparam;
  std::vector<float> values(nvalues);
  deSerialize<float>(stream, nvalues, &values[0]);

  cfgScalar scalingFactor = 290.91/200;

  //Relative Young's modulus,	Poisson's ratio,	Relative shear modulus,	Density
  assert(values.size()%dimparam==0);
  int ival, nval=(int)values.size();
  for (ival=0; ival<nval; ival+=dimparam)
  {
    cfgScalar E= scalingFactor*values[ival+1];
    cfgScalar nu = values[ival+2];
    cfgScalar mu = E/(2*(1+nu));
    cfgScalar density = 1.;

    ioPhysicalParameters.push_back(density);
    ioPhysicalParameters.push_back(E);
    ioPhysicalParameters.push_back(nu);
    ioPhysicalParameters.push_back(mu);
  }

  return true;
}

std::vector<std::vector<int> > SamplesGeneratorImpl::getBaseMaterialStructures(int nmat, int icell)
{
  assert(icell==1);

  std::vector<std::vector<int> > baseMaterialStructures(nmat);
  for (int imat=0; imat<nmat; imat++)
  {
    baseMaterialStructures[imat].push_back(imat);
  }
  return baseMaterialStructures;
}

int SamplesGeneratorImpl::runExhaustiveGamutComputation(int iLevel)
{
  assert(iLevel <=2 || ( iLevel<=8 && (m_cubicOnly || m_orthotropicOnly) ) );

  int level = iLevel;
  int n[3] = {level, level, level};

  std::string stepperType = "newton";
  std::vector<std::vector<int> > newMaterialAssignment;
  computeMaterialParametersIncremental(stepperType, level, newMaterialAssignment);

  std::vector<float> tensors;
  std::vector<float> physicalParameters;
  computeParametersAndTensorValues(n, newMaterialAssignment, physicalParameters, tensors);
  writeFiles(level, newMaterialAssignment, getBaseMaterialStructures(m_nmat), physicalParameters, tensors);

  return 0;
}

int SamplesGeneratorImpl::runDiscreteAndContinuousOptimization(int iLevel, int iNbCycles, int iStartCycle, bool iGrowStructure)
{
  int level = iLevel;
  int prevLevel = level/2;
  int growStructure = iGrowStructure;
  bool fixNonManifoldStructure = true;
  std::string stepperType = "newton";

  int icycle=0, ncycle=iNbCycles;
  for (icycle=1; icycle<ncycle; icycle++)
  {
    if (icycle>0)
    {
      growStructure = false;
    }
    std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
    std::vector<float> physicalParameters, tensors;
    bool ResOk = readFiles((growStructure? prevLevel: level), materialAssignments, baseMaterialStructures, physicalParameters, tensors);
    if (ResOk)
    {
      std::vector<std::vector<int> > allNewMaterialAssignments;
      std::vector<cfgScalar> allParameters, allTensors;
      std::vector<std::vector<std::vector<float> > > x[2];

      computeVariationsSMC(level, stepperType, materialAssignments, baseMaterialStructures, physicalParameters, tensors, growStructure, allNewMaterialAssignments, allParameters, allTensors); 
      writeFiles(level, allNewMaterialAssignments, baseMaterialStructures, allParameters, allTensors, "SMC_" +std::to_string(icycle));
      //writeFiles(level, allNewMaterialAssignments, baseMaterialStructures, allParameters, allTensors, "_" +std::to_string(icycle));

      std::vector<std::vector<int> > newMaterialAssignments;
      std::vector<cfgScalar> newParameters, newTensors;
      runContinuousOptimization(level, icycle, fixNonManifoldStructure, allNewMaterialAssignments, allParameters, allTensors, newMaterialAssignments, newParameters, newTensors);
      writeFiles(level, newMaterialAssignments, baseMaterialStructures, newParameters, newTensors, "ContinuousOptimizationResult_"+std::to_string(icycle));

      writeFiles(level, newMaterialAssignments, baseMaterialStructures, newParameters, newTensors);
    }
  }
  return 0;
}

int SamplesGeneratorImpl::runDiscreteOptimization(int iLevel, int iNbCycles, int iStartCycle, bool iGrowStructure)
{
  int level = iLevel;
  int prevLevel = level/2;
  bool growStructure = iGrowStructure;
  std::string stepperType = "newton";

  int icycle=0, ncycle=iNbCycles;
  for (icycle=iStartCycle; icycle<ncycle; icycle++)
  {
    std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
    std::vector<float> physicalParameters, tensors;
    std::string fileName = "";
    if (!growStructure && icycle>0)
    {
      fileName = "SMC_" +std::to_string(icycle-1);
    }
    bool ResOk = readFiles((growStructure? prevLevel: level), materialAssignments, baseMaterialStructures, physicalParameters, tensors, fileName);
    if (ResOk)
    {
      std::vector<std::vector<int> > allNewMaterialAssignments;
      std::vector<cfgScalar> allParameters, allTensors;
      std::vector<std::vector<std::vector<float> > > x[2];

      computeVariationsSMC(level, stepperType, materialAssignments, baseMaterialStructures, physicalParameters, tensors, growStructure, allNewMaterialAssignments, allParameters, allTensors); 
      writeFiles(level, allNewMaterialAssignments, baseMaterialStructures, allParameters, allTensors, "SMC_" +std::to_string(icycle));

      writeFiles(level, allNewMaterialAssignments, baseMaterialStructures, allParameters, allTensors);
    }
    growStructure = false;
  }
  return 0;
}

int SamplesGeneratorImpl::runStrengthAnalysis(int iLevel)
{
  int level = iLevel;
  std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
  std::vector<float> physicalParameters, tensors;
  bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters, tensors);
  if (ResOk)
  {
    int n[3] = {level, level, m_dim==2?1: level};
    StrengthAnalysis::StructureType type = (m_cubicOnly? StrengthAnalysis::Cubic: (m_orthotropicOnly? StrengthAnalysis::Orthotropic:  StrengthAnalysis::General));

    StrengthAnalysis strengthAnalysis;
    std::vector<float> strengths;
    if (m_dim==2)
    {
      strengthAnalysis.run2D(n, type, m_mat2D, materialAssignments, physicalParameters, tensors, m_blockRep, m_nbSubdivisions, strengths);
    }
    else
    {
      strengthAnalysis.run3D(n, type, m_mat, materialAssignments, physicalParameters, tensors, m_blockRep, m_nbSubdivisions, strengths);
    }
    if (strengths.size()>0)
    {
    std::string fileRootName = m_OutputDirectory + "level" + std::to_string(level) + "_";
    std::string fileExtension = ".bin";
    cfgUtil::writeBinary<float>(fileRootName + "ultimateStrengths" + fileExtension, strengths);
    }
  }
  return 0;
}


int SamplesGeneratorImpl::runVariantsComputation(int iLevel)
{
  int level = iLevel;
  std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
  std::vector<float> physicalParameters, tensors;
  bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters, tensors);
  if (ResOk)
  {
    std::vector<std::vector<int> > newMaterialAssignments;
    int n[3] = {level, level, m_dim==2?1: level};
    int t[3] = {n[0]/2, n[1]/2, n[2]/2};
    std::set<std::vector<int> > setMaterialAssignments;
    std::vector<int> indices;
    //std::map<std::vector<int>, int> map;

    int imat, nmat=(int)materialAssignments.size();
    for (imat=0; imat<nmat; imat++)
    {
      const std::vector<int> & matAssignment = materialAssignments[imat];
      if (setMaterialAssignments.count(matAssignment)==0)
      {
        newMaterialAssignments.push_back(matAssignment);
        indices.push_back(imat);
        setMaterialAssignments.insert(matAssignment);
        //map[matAssignment] = imat;

        std::vector<int> translatedStructure;
        if (m_dim==2)
        {
          translateStructure(n[0], n[1], matAssignment, t[0], t[1], translatedStructure);
        }
        else
        {
          translateStructure(n[0], n[1], n[2], matAssignment, t[0], t[1], t[2], translatedStructure);
        }
        if (setMaterialAssignments.count(translatedStructure)==0)
        {
          newMaterialAssignments.push_back(translatedStructure);
          indices.push_back(imat);
          setMaterialAssignments.insert(translatedStructure);
        } 
      }
      /*else
      {
        std::vector<int> ind1, ind2;
        ind1.push_back(imat);
        int oldMatIndex = map[matAssignment];
        ind2.push_back(oldMatIndex);
        std::vector<float> initParams, newParams;
        initParams = getSubVector(physicalParameters, getNbParameters(), ind1);
        newParams = getSubVector(physicalParameters, getNbParameters(), ind2);
      }*/ 
    }
    std::vector<float> newPhysicalParameters = getSubVector(physicalParameters, getNbParameters(), indices);
    std::vector<float> newTensors = getSubVector(tensors, (m_dim==2?6:21), indices);
    writeFiles(level, newMaterialAssignments, baseMaterialStructures, newPhysicalParameters, newTensors, "with_variants");   
  }
  return 0;
}

int SamplesGeneratorImpl::runFilesConcatenation()
{
  int level = 16;
  int indexMin = 0;
  int indexMax = 4;
  //std::string postFix = "SMC_after_resampling_";
  //std::string newPostFix = "SMC_result";
  std::string postFix = "ContinuousOptimizationResult_";
  std::string newPostFix = "ContinuousOptimizationResult";
  //std::string postFix = "SMC_";
  //std::string newPostFix = "SMC";
  //concatenateFiles(level, indexMin, indexMax, postFix, newPostFix);

  std::string postFix1 = "";
  std::string postFix2 = "3D_fixed";
  //std::string postFix1 = "SMC_init";
  //std::string postFix2 = "SMC_result";
  //std::string postFix1 = "ContinuousOptimizationResult";
  //std::string postFix2 = "SMC";
  newPostFix = "";
  concatenateFiles(level, postFix1, postFix2, newPostFix);
  return 0;
}

int SamplesGeneratorImpl::runIsotropicToOrthotropicConversion(int iLevel)
{
  int level = iLevel;
  std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
  std::vector<float> cubicParameters, tensors;
  bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, cubicParameters, tensors);

  std::vector<float> orthotropicParameter;
  fromCubicToOrthotropicProperties(m_dim, cubicParameters, orthotropicParameter);
  writeFiles(level, materialAssignments, baseMaterialStructures, orthotropicParameter, tensors, "orthotropic");
  return 0;
}

int SamplesGeneratorImpl::runSubsamplingStage(int iLevel, const std::string &iSuffix)
{
  int level = iLevel;
  std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
  std::vector<float> physicalParameters, physicalParametersInit, tensors;
  std::string fileName = iSuffix;
  bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters, tensors, fileName);
  if (!ResOk)
    return 1;

  physicalParametersInit = physicalParameters;

  int paramdim = getNbParameters();
  if (m_useLogScale)
  {
    std::vector<int> dimToScale;
    if (1)
    {
      if (paramdim==4)
      {
        dimToScale.push_back(1);
      }
      else if (paramdim==5)
      {
        dimToScale.push_back(1);
        dimToScale.push_back(2);
      }
      else
      {
        dimToScale.push_back(1);
        dimToScale.push_back(2);
        dimToScale.push_back(3);
      }
    }
    std::vector<cfgScalar> lengths(paramdim, 1);
    rescaleData(physicalParameters, paramdim, lengths);

    cfgScalar eps = 1.e-6;
    convertToLogValues(physicalParameters, paramdim, dimToScale, eps);
  }

  int nTargetParticules = 50000;
  //cfgScalar minRadius = 0.02;
  cfgScalar minRadius = 0.2;

  std::vector<int> particules;

  bool greedy = false;
  if (greedy)
  {
    //samplePointsRandom(nTargetParticules, physicalParameters, paramdim, minRadius, particules);
    samplePointsGreedyV2(nTargetParticules, physicalParameters, paramdim, minRadius, particules);
    //getFurthestPointsGreedy(nTargetParticules, physicalParameters, paramdim, minRadius, particules);
  }
  else
  {
    bool boundaryOnly = false;

    if (boundaryOnly)
    {
      DistanceField distanceField(paramdim);
      std::vector<cfgScalar> distances = distanceField.computeDistances(physicalParameters);

      Resampler resampler;
      resampler.resampleBoundary(minRadius, paramdim, physicalParameters, distances, nTargetParticules, particules);
      std::cout << "nb points = " << particules.size() << std::endl;
    }
    else
    {
      if (0)
      {
        ScoringFunction scoringFunction(paramdim);
        bool useDistanceField = true;
        scoringFunction.setUseDistanceField(useDistanceField);
        std::vector<cfgScalar> scores = scoringFunction.computeScores(physicalParameters);

        Resampler resampler;
        resampler.resample(minRadius, paramdim, physicalParameters, scores, nTargetParticules, particules);
      }
      else
      {
        DistanceField distanceField(paramdim);
        std::vector<cfgScalar> distances = distanceField.computeDistances(physicalParameters);

        Resampler resampler;
        resampler.resampleUsingNormalDistribution(distances, nTargetParticules, particules);
      }
    }
  }
  std::vector<cfgScalar> initParameters, initTensors;
  physicalParametersInit = getSubVector(physicalParametersInit, paramdim, particules);
  tensors = getSubVector(tensors,  (m_dim==2?6:21), particules);
  materialAssignments = getSubVector(materialAssignments, particules);
  writeFiles(level, materialAssignments, baseMaterialStructures, physicalParametersInit, tensors, iSuffix+"_subsampled");

  return 0;
}

int SamplesGeneratorImpl::runMirroring(int iLevel)
{

  // mirror structures
  int level = iLevel;
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
  return 0;
}

int SamplesGeneratorImpl::run2Dto3DConversion(int iLevel)
{
  // Convert 2D cubic structures to 3D cubic structures
  int level = iLevel;
  bool fixNonManifoldStructures = true;
  int n[3] = {level, level, level};
  std::vector<std::vector<int> > materialAssignments2D, baseMaterialStructures;
  std::vector<float> cubicParameters, tensors;
  bool ResOk = readFiles(level, materialAssignments2D, baseMaterialStructures, cubicParameters, tensors, "subsampled");

  std::vector<std::vector<int> > materialAssignments3D;
  int imat, nmat=(int)materialAssignments2D.size();
  materialAssignments3D.resize(nmat);
  for (imat=0; imat<nmat; imat++)
  {
    convert2DCubicStructuresTo3DCubicStructures(level, materialAssignments2D[imat], materialAssignments3D[imat]);
    if (fixNonManifoldStructures)
    {
      fixNonManifoldStructure(n, materialAssignments3D[imat]);
    }
  }
  std::vector<float> cubicParameters3D, tensors3D;
  computeParametersAndTensorValues(n, materialAssignments3D, cubicParameters3D, tensors3D);
  writeFiles(level, materialAssignments3D, baseMaterialStructures, cubicParameters3D, tensors3D, "3D");
  return 0;
}

int SamplesGeneratorImpl::runUpscaleStructure(int iLevel)
{
  int level = iLevel;
  bool recomputeParam = true;

  std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
  std::vector<float> physicalParameters, tensors;
  bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters, tensors);

  int newlevel = 2*level;
  std::vector<std::vector<int> > newMaterialAssignments;
  int imat, nmat=(int)materialAssignments.size();
  for (imat=0; imat<nmat; imat++)
  {
    std::vector<int> newMatAssignment;
    if (m_dim==2)
    {
      upscaleStructure(level, level, materialAssignments[imat], newMatAssignment);
    }
    else
    {
      upscaleStructure(level, level, level, materialAssignments[imat], newMatAssignment);
    }
    newMaterialAssignments.push_back(newMatAssignment);
  } 
  if (recomputeParam)
  {
    int n[3] = {newlevel, newlevel, newlevel};
    std::vector<float> newtensors;
    std::vector<float> newphysicalParameters;
    computeParametersAndTensorValues(n, newMaterialAssignments, newphysicalParameters, newtensors);
    writeFiles(newlevel, newMaterialAssignments, baseMaterialStructures, newphysicalParameters, newtensors);
  }
  else
  {
    writeFiles(newlevel, newMaterialAssignments, baseMaterialStructures, physicalParameters, tensors);
  }
  return 0;
}

int SamplesGeneratorImpl::runMaterialPropertiesComputation(int iLevel, const std::string &iSuffix)
{
  int level = iLevel;
  int n[3] = {level, level, level};

  std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
  std::vector<float> physicalParameters, tensors;
  bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters, tensors, iSuffix);

  std::string suffix = iSuffix+"_updated";
  computeParametersAndTensorValues(n, materialAssignments, physicalParameters, tensors);
  writeFiles(level, materialAssignments, getBaseMaterialStructures(m_nmat), physicalParameters, tensors, suffix);

  return 0;
}

int SamplesGeneratorImpl::runTensorHomogenization(int iLevel)
{
  int level = iLevel;
  std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
  std::vector<float> physicalParameters;
  bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters);
  
  //std::vector<std::vector<double> > materialDistributions;
  //bool ResOk = cfgUtil::readBinary<double>(m_OutputDirectory + "16Ortho.bin", materialDistributions);

  int n[3] = {level, level, level};
  std::vector<float> tensors;
  std::vector<float> newphysicalParameters;
  computeParametersAndTensorValues(n, materialAssignments, newphysicalParameters, tensors);
  writeFiles(level, materialAssignments, baseMaterialStructures, newphysicalParameters, tensors);
  //writeFiles(level, convertVec<double, cfgScalar>(materialDistributions), physicalParameters, tensors, "continuous");

  return 0;
}

int SamplesGeneratorImpl::runFamilyGeneration(int iLevel)
{
  int level = iLevel;
  int n[3] = {level, level, level};

  FamilyGenerator * generator = FamilyGenerator::createOperator();
  generator->setMicrostructureSize(m_dim, n);
  generator->run();
  const std::vector<std::vector<int> > & matAssignments = generator->getMicrostructures();

  std::vector<float> physicalParameters;
  std::string suffix = "newFamily";
  std::vector<float> tensors;
  computeParametersAndTensorValues(n, matAssignments, physicalParameters, tensors);
  writeFiles(level, matAssignments, getBaseMaterialStructures(m_nmat), physicalParameters, tensors, suffix);

  delete generator;

  return 0;
}

// thermal properties
int SamplesGeneratorImpl::runThermalAnalysis(int iLevel)
{
  int level = iLevel;
  std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
  std::vector<float> physicalParameters, tensors;
  bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters, tensors);
  if (ResOk)
  {
    int n[3] = {level, level, m_dim==2?1: level};

    ThermalAnalysis thermalAnalysis;
    std::vector<cfgScalar> thermalCoefficients;
    if (m_dim==2)
    {
      thermalAnalysis.run2D(n, m_mat2D, materialAssignments, physicalParameters, tensors, m_blockRep, m_nbSubdivisions, thermalCoefficients);
    }
    else
    {
      thermalAnalysis.run3D(n, m_mat, materialAssignments, physicalParameters, tensors, m_blockRep, m_nbSubdivisions, thermalCoefficients);
    }

    std::string fileRootName = m_OutputDirectory + "level" + std::to_string(level) + "_";
    std::string fileExtension = ".bin";
    cfgUtil::writeBinary<float>(fileRootName + "thermalExpansionCoefficients" + fileExtension, thermalCoefficients);
  }
  return 0;
}

int SamplesGeneratorImpl::runFixDisconnectedMicrostructures(int iLevel, const std::string &iSuffix)
{
  int level = iLevel;
  int n[3] = {level, level, level};
  std::vector<int> connectedStructures, fixedStructures, invalidStructures;
  std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
  std::vector<float> physicalParameters, tensors;

  bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters, tensors, iSuffix);
  if (tensors.size()==0)
  {
    int nmat=(int)materialAssignments.size();
    int tensorDim = (m_dim==2?6:21);
    tensors.resize(tensorDim*nmat);
  }

  if (ResOk)
  {
    std::vector<int> disconnectedStructures;
    int imat, nmat=(int)materialAssignments.size();
    std::cout << "nmat init = " << nmat << std::endl;
    for (imat=0; imat<nmat; imat++)
    {
      if (m_dim==2)
      {
        if (imat % 10000 == 1)
        {
          std::cout << "imat = " << imat << std::endl;
        }
      }
      else
      {
        if (imat % 1000 == 1)
        {
          std::cout << "imat = " << imat << std::endl;
        }
      }

      bool modified = false;
      bool resOk = true;
      if (m_dim==2)
      {
        resOk = filterOutNonConnectedComponents(n[0], n[1], materialAssignments[imat], modified);
      }
      else
      {
        resOk = filterOutNonConnectedComponents(n[0], n[1], n[2], materialAssignments[imat], modified);
      }
      if (resOk)
      {
        if (modified)
        {
          fixedStructures.push_back(imat);
        }
        else
        {
          connectedStructures.push_back(imat);
        }
      }
      else
      {
        invalidStructures.push_back(imat);
      }
    }
    std::vector<std::vector<int> > fixedMaterialAssignments = getSubVector(materialAssignments, fixedStructures);
    std::vector<float> fixedPhysicalParameters = getSubVector(physicalParameters, getNbParameters(), fixedStructures);
    std::vector<float> fixedTensors = getSubVector(tensors, (m_dim==2?6:21), fixedStructures);
    if (fixedMaterialAssignments.size()>0)
    {
      std::cout << "nb fixed_structures = " << fixedMaterialAssignments.size() << std::endl;
      writeFiles(level, fixedMaterialAssignments, baseMaterialStructures, fixedPhysicalParameters, fixedTensors, "fixed_structures");
    }
    std::vector<std::vector<int> > invalidMaterialAssignments = getSubVector(materialAssignments, invalidStructures);
    std::vector<float> invalidPhysicalParameters = getSubVector(physicalParameters, getNbParameters(), invalidStructures);
    std::vector<float> invalidTensors = getSubVector(tensors, (m_dim==2?6:21), invalidStructures);
    if (invalidMaterialAssignments.size()>0)
    {
      std::cout << "nb invalid_structures = " << invalidMaterialAssignments.size() << std::endl;
      writeFiles(level, invalidMaterialAssignments, baseMaterialStructures, invalidPhysicalParameters, invalidTensors, "invalid_structures");
    }
    std::vector<std::vector<int> > validMaterialAssignments = getSubVector(materialAssignments, connectedStructures);
    std::vector<float> validPhysicalParameters = getSubVector(physicalParameters, getNbParameters(), connectedStructures);
    std::vector<float> validTensors = getSubVector(tensors, (m_dim==2?6:21), connectedStructures);
    if (validMaterialAssignments.size()>0)
    {
      std::cout << "nb connected_structures = " << validMaterialAssignments.size() << std::endl;
      writeFiles(level, validMaterialAssignments, baseMaterialStructures, validPhysicalParameters, validTensors, "connected_structures");
    }
    std::vector<float> newtensors;
    std::vector<float> newphysicalParameters;
    std::cout << "fixed mat = " << fixedMaterialAssignments.size() << std::endl;
    if (fixedMaterialAssignments.size()>0)
    {
      computeParametersAndTensorValues(n, fixedMaterialAssignments, newphysicalParameters, newtensors);

      validMaterialAssignments.insert(validMaterialAssignments.end(), fixedMaterialAssignments.begin(), fixedMaterialAssignments.end());
      validPhysicalParameters.insert(validPhysicalParameters.end(), newphysicalParameters.begin(), newphysicalParameters.end());
      validTensors.insert(validTensors.end(), newtensors.begin(), newtensors.end());
    }
    if (validMaterialAssignments.size()>0)
      writeFiles(level, validMaterialAssignments, baseMaterialStructures, validPhysicalParameters, validTensors, "valid");
  }
  return 0;
}

int SamplesGeneratorImpl::runFixNonManifoldMicrostructure(int iLevel)
{
  int level = iLevel;
  int n[3] = {level, level, level};

  std::vector<std::vector<int> > materialAssignments, baseMaterialStructures;
  std::vector<float> physicalParameters, tensors;
  bool ResOk = readFiles(level, materialAssignments, baseMaterialStructures, physicalParameters, tensors, "3D");
  if (ResOk)
  {
    std::vector<int> nonManifoldStructuresIndices;
    int imat, nmat=(int)materialAssignments.size();
    for (imat=0; imat<nmat; imat++)
    {
      bool fixed = fixNonManifoldStructure(n, materialAssignments[imat]);
      if (fixed)
      {
        nonManifoldStructuresIndices.push_back(imat);
      }
    }
    std::vector<std::vector<int> > nonManifoldMaterialAssignments = getSubVector(materialAssignments, nonManifoldStructuresIndices);
    std::vector<float> nonManifoldPhysicalParameters = getSubVector(physicalParameters, getNbParameters(), nonManifoldStructuresIndices);
    std::vector<float> nonManifoldTensors = getSubVector(tensors, (m_dim==2?6:21), nonManifoldStructuresIndices);
    writeFiles(level, nonManifoldMaterialAssignments, baseMaterialStructures, nonManifoldPhysicalParameters, nonManifoldTensors, "3D_nonmanifold");

    std::vector<float> newtensors;
    std::vector<float> newphysicalParameters;
    computeParametersAndTensorValues(n, nonManifoldMaterialAssignments, newphysicalParameters, newtensors);

    setSubVector(physicalParameters, newphysicalParameters, getNbParameters(), nonManifoldStructuresIndices);
    setSubVector(tensors, newtensors, (m_dim==2?6:21), nonManifoldStructuresIndices);
    writeFiles(level, materialAssignments, baseMaterialStructures, physicalParameters, tensors, "3D_fixed");
  }
  return 0;
}

int SamplesGeneratorImpl::rewriteDisneyFiles()
{
  int dim = 2;
  int cubic = false;
  int orthotropic = true;
  std::vector<cfgScalar> physicalParameters;
  bool resOk = readDisneyFiles(dim, cubic, orthotropic, physicalParameters);
  if (resOk)
  {
    std::string fileRootName = m_OutputDirectory + "level" + std::to_string(0) + "_" + "disney_2D_orthotropic_";
    std::string fileExtension = ".bin";
    cfgUtil::writeBinary<float>(fileRootName + "params" + fileExtension, physicalParameters);
  }
  return 0;
}

int SamplesGeneratorImpl::rewritePanettaFiles()
{
  int dim = 3;
  std::vector<cfgScalar> physicalParameters;
  bool resOk = readPanettaFiles(dim, physicalParameters);
  if (resOk)
  {
    std::string fileRootName = m_OutputDirectory + "level" + std::to_string(0) + "_" + "Panetta_3D_cubic_with_rescaled_E_";
    std::string fileExtension = ".bin";
    cfgUtil::writeBinary<float>(fileRootName + "params" + fileExtension, physicalParameters);
  }
  return 0;
}

void SamplesGeneratorImpl::init()
{
  std::vector<float> E(m_nmat), nu(m_nmat), lambda(m_nmat), mu(m_nmat);

  bool setBaseStiffness = false;
  bool updatePoissonRatio = true;

  if (setBaseStiffness)
  {
    E[1] = 200;
    E[0] = E[1]/1000;

    nu[1] = 0.35;
    nu[0] = 0.35;

    fromYoungModulusPoissonRatioToLamesParameters(E[0], nu[0], lambda[0], mu[0]);
    fromYoungModulusPoissonRatioToLamesParameters(E[1], nu[1], lambda[1], mu[1]);
  }
  else
  {
    mu[1] = 100;
    lambda[1] = 1000;
    fromLamesParametersToYoungModulusPoissonRatio(lambda[1], mu[1], E[1], nu[1]);

    E[0] = E[1]/1000;
    if (m_nmat==2)
    {
      nu[0] = 0.1;
    }
    else
    {
      nu[0] = nu[1];
    }
    fromYoungModulusPoissonRatioToLamesParameters(E[0], nu[0], lambda[0], mu[0]);
    if (updatePoissonRatio)
    {
      nu[0] = 0.48;
      nu[1] = nu[0];

      fromYoungModulusPoissonRatioToLamesParameters(E[0], nu[0], lambda[0], mu[0]);
      fromYoungModulusPoissonRatioToLamesParameters(E[1], nu[1], lambda[1], mu[1]);
    }
  }
  if (m_nmat>2)
  {
    mu[2] = 0;
    lambda[2] = 0;
  }
  m_baseMaterialDensities.clear();
  m_baseMaterialDensities.resize(m_nmat, 1);
  m_baseMaterialDensities[m_nmat-1] = 0;

  if (m_nmat==2)
  {
    //ratio of rigid
    std::swap(m_baseMaterialDensities[0], m_baseMaterialDensities[1]);
  }

  m_ene.resize(m_nmat);
  m_ene2D.resize(m_nmat);
  for (int imat=0; imat<m_nmat; imat++)
  {
    m_ene2D[imat].param[0] = mu[imat];
    m_ene2D[imat].param[1] = lambda[imat];

    m_ene[imat].param[0] = mu[imat];
    m_ene[imat].param[1] = lambda[imat];
  }

  if (m_dim==2)
  {
    m_mat2D.clear();
    m_mat2D.resize(m_ene2D.size());
    for (unsigned int ii = 0; ii < m_mat2D.size(); ii++){
      for (unsigned int jj = 0; jj < m_mat2D[ii].e.size(); jj++){
        m_mat2D[ii].e[jj] = &m_ene2D[ii];
      }
    }
  }
  else
  { 
    m_mat.clear();
    m_mat.resize(m_ene.size());
    for (unsigned int ii = 0; ii < m_mat.size(); ii++){
      for (unsigned int jj = 0; jj < m_mat[ii].e.size(); jj++){
        m_mat[ii].e[jj] = &m_ene[ii];
      }
    }
  }

  m_parametersToUse.clear();
  bool useAllParameters = false;
  if (useAllParameters)
  {
    int nparam = getNbParameters();
    m_parametersToUse = genIncrementalSequence(0, nparam-1);
  }
  else
  {
    if (m_dim==2)
    {
      m_parametersToUse.push_back(0); // density
      if (m_cubicOnly)
      {
        m_parametersToUse.push_back(1); // E
        m_parametersToUse.push_back(2); //nu
      }
      else if (m_orthotropicOnly)
      {
        m_parametersToUse.push_back(1); // Ex
        m_parametersToUse.push_back(2); //Ey
        m_parametersToUse.push_back(3); //Nu_xy
      }
      else
      {
        assert(0);
      }
    }
    else
    {
      if (m_cubicOnly)
      {
        m_parametersToUse.push_back(1); // E
        m_parametersToUse.push_back(2); //nu
      }
      else
      {
        assert(0);
      }
    }
  }
}

int SamplesGeneratorImpl::run()
{
  // Init parameters
  Stage stage = DiscreteOptimizationStage;
  int level = 16;
  int ncycle = 1;
  int startCycle = 0;
  bool growStructure = false;

  std::string suffix = "";

  const char * filename = "config.txt";
  ConfigFile conf;
  conf.load(filename);

  if (conf.hasOpt("suffix"))
  {
    std::string str = conf.getString("suffix");
    suffix = str;
  }
  if (conf.hasOpt("type"))
  {
    std::string str = conf.getString("type");
    if (str == "Cubic2D")
    {
      m_dim = 2;
      m_cubicOnly = true;
      setOutputDirectory("..//..//Output_2D_cubic//");
    }
    else if (str == "Orthotropic2D")
    {
      m_dim = 2;
      m_orthotropicOnly = true;
      setOutputDirectory("..//..//Output_2D_ortho//");
    }
    else if (str == "Cubic3D")
    {
      m_dim = 3;
      m_cubicOnly = true;
      setOutputDirectory("..//..//Output_3D_cubic//");
    }
    else
    {
      assert(0);
    }
  }
  if (conf.hasOpt("outputDir"))
  {
    std::string outputDir = conf.getString("outputDir");
    setOutputDirectory(outputDir);
  }
  if (conf.hasOpt("nmat"))
  {
    m_nmat = conf.getInt("nmat");;
  }
  if (conf.hasOpt("stage"))
  {
    std::string str = conf.getString("stage");
    if (str == "ExhaustiveGamutComputation")
    {
      stage = ExhaustiveGamutComputationStage;
    }
    else if (str == "DiscreteAndContinuousOptimization")
    {
      stage = DiscreteAndContinuousOptimizationStage;
    }
    else if (str == "DiscreteOptimization")
    {
      stage = DiscreteOptimizationStage;
    }
    else if (str == "ContinuousOptimization")
    {
      stage = ContinuousOptimizationStage;
    }
    else if (str == "StrengthAnalysis")
    {
      stage = StrengthAnalysisStage;
    }
    else if (str == "VariantsComputation")
    {
      stage = VariantsComputationStage;
    }
    else if (str == "FilesConcatenation")
    {
      stage = FilesConcatenationStage;
    }
    else if (str == "IsotropicToOrthotropicConversion")
    {
      stage = IsotropicToOrthotropicConversionStage;
    }
    else if (str == "Subsampling")
    {
      stage = SubsamplingStage;
    }
    else if (str == "MirrorStructure")
    {
      stage = MirrorStructureStage;
    }
    else if (str == "2Dto3DConversion")
    {
      stage = From2Dto3DConversionStage;
    }
    else if (str == "UpscaleMicrostructure")
    {
       stage = UpscaleMicrostructureStage;
    }
    else if (str == "FixDisconnectedMicrostructure")
    {
      stage = FixDisconnectedMicrostructureStage;
    }
    else if (str == "FixNonManifoldMicrostructure")
    {
      stage = FixNonManifoldMicrostructureStage;
    }
    else if (str == "MaterialPropertiesComputation")
    {
      stage = MaterialPropertiesComputationStage;
    }
    else if (str == "homogenizeTensor")
    {
      stage = TensorHomogenizationStage;
    }
    else if (str == "FamilyGeneration")
    {
      stage = FamilyGenerationStage;
    }
    else
    {
      std::cout << "invalid stage!" << std::endl;
      assert(0);
    }
  }
  if (conf.hasOpt("level"))
  {
    level = conf.getInt("level");
  }
  if (conf.hasOpt("nbCycles"))
  {
    ncycle = conf.getInt("nbCycles");
  }
  if (conf.hasOpt("startCycle"))
  {
    startCycle = conf.getInt("startCycle");
  }
  if (conf.hasOpt("growStructure"))
  {
    growStructure = conf.getBool("growStructure");
  }
  if (conf.hasOpt("filterOutDisconnetedStructures"))
  {
    m_filterOutDisconnetedStructures = conf.getBool("filterOutDisconnetedStructures");
  }
  if (conf.hasOpt("useLogScale"))
  {
    m_useLogScale = conf.getBool("useLogScale");
  }

  init();

  srand(0);

  if (stage == ExhaustiveGamutComputationStage)
  {
    int status = runExhaustiveGamutComputation(level);
    return status;
  }
  else if (stage == DiscreteAndContinuousOptimizationStage)
  {
    int status = runDiscreteAndContinuousOptimization(level, ncycle, startCycle, growStructure);
    return status;
  }
  else if (stage == DiscreteOptimizationStage) // SMC
  {
    int status = runDiscreteOptimization(level, ncycle, startCycle, growStructure);
    return status;
  }
  else if (stage == ContinuousOptimizationStage)
  {
    int status = runContinuousOptimization(level, startCycle, suffix);
    return status;
  }
  else if (stage == StrengthAnalysisStage)
  {
    int status = runStrengthAnalysis(level);
    return status;
  }
  else if (stage == VariantsComputationStage)
  {
    int status = runVariantsComputation(level);
    return status;
  }
  else if (stage == FilesConcatenationStage)
  {
    int status = runFilesConcatenation();
    return status;
  }
  else if (stage == IsotropicToOrthotropicConversionStage)
  {
    int status = runIsotropicToOrthotropicConversion(level);
    return status;
  }
  else if (stage == SubsamplingStage)
  {
    int status = runSubsamplingStage(level, suffix);
    return status;
  }
  else if (stage == MirrorStructureStage)
  {
    int status = runMirroring(level);
    return status;
  }
  else if (stage == From2Dto3DConversionStage)
  {
    int status = run2Dto3DConversion(level);
    return status;
  }
  else if (stage == UpscaleMicrostructureStage)
  {
    int status = runUpscaleStructure(level);
  }
  else if (stage == FixDisconnectedMicrostructureStage)
  {
    int status = runFixDisconnectedMicrostructures(level, suffix);
    return status;
  }
  else if (stage == FixNonManifoldMicrostructureStage)
  {
    int status = runFixNonManifoldMicrostructure(level);
  }
  else if (stage == MaterialPropertiesComputationStage)
  {
    int status = runMaterialPropertiesComputation(level, suffix);
  }
  else if (stage == TensorHomogenizationStage)
  {
    int status = runTensorHomogenization(level);
  }
  else if (stage == FamilyGenerationStage)
  {
    int status = runFamilyGeneration(level);
  }
  else if (stage == ThermalAnalysisStage)
  {
    int status = runThermalAnalysis(level);
  }
  return 0;
}




