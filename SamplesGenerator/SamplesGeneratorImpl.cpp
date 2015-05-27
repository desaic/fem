
#include "SamplesGeneratorImpl.h"

#include <assert.h>

#include "ElementRegGrid.hpp"
#include "Element.hpp"

#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"

#include "Stepper.hpp"
#include "StepperNewton.hpp"
#include "AdmmCPU.hpp"
#include "StepperGrad.hpp"

#include "cfgUtilities.h"
using namespace cfgUtil;

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;


SamplesGeneratorImpl::SamplesGeneratorImpl()
{
}

SamplesGeneratorImpl::~SamplesGeneratorImpl()
{
}

float SamplesGeneratorImpl::computeStrain(const ElementRegGrid * iElementGrid, int iAxis)
{
  writeVector2File(iElementGrid->x, "Output//x.m");

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
}

bool SamplesGeneratorImpl::computeEquilibrium(Stepper * iStepper, ElementRegGrid * iElementGrid, int iNumberOfSteps)
{
  assert(iStepper);

  iStepper->nSteps = iNumberOfSteps;
  iStepper->init(iElementGrid);
  int status = iStepper->step();
  bool Converged = status<0;
  return Converged;
}

void SamplesGeneratorImpl::computeForceDeformationSample(ElementRegGrid * iElementGrid, std::string iStepperType, int iNumberOfSteps)
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
  bool Converged = computeEquilibrium(stepper, iElementGrid, iNumberOfSteps);
  
  delete(stepper);
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

void SamplesGeneratorImpl::sampleMaterialSpace()
{
  std::string fileDirectory = "C://Users//melina//Documents//Projects//fem//vc//runtime//Output//";

  int level = 2;
  int blockSize = level+1;

  std::string materialFile = fileDirectory + "Materials_" + std::to_string(level) + ".txt";
  std::vector<std::vector<int> > materialCombinations;
  readMaterialCombinations(materialFile, materialCombinations);

  std::vector<Scalar> physicalParameters; //YoungModulus, density;
  std::vector<std::vector<int> > macroCellMaterials;
  int icomb, ncomb=256;
  macroCellMaterials.resize(ncomb);
  for (icomb=0; icomb<ncomb; icomb++)
  {
    std::string sampleFile = fileDirectory + "StressStrain_" + std::to_string(level) + "_" + std::to_string(blockSize) + "_" + std::to_string(icomb) + ".txt"; 
    Vector3f forceAxis;
    std::vector<int> materials;
    std::vector<Scalar> stresses, strains;
    readData(sampleFile, forceAxis, materials, stresses, strains);

    Scalar r=0;
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
    r /= (Scalar)ncell;

    Scalar Y = computeYoungModulus(strains, stresses);
    physicalParameters.push_back(Y);
    physicalParameters.push_back(r);
  }
  std::vector<int> convexHullVertices;
  computeConvexHull(physicalParameters, 2, convexHullVertices);
  std::vector<std::vector<int> > newMaterialCombinations = getSubVector(macroCellMaterials, convexHullVertices);
}

void SamplesGeneratorImpl::sampleDeformation(int iN[3], std::vector<MaterialQuad> &iMaterial, const std::string iStepperType, float iMaxForce, int iForceAxis, int iNumberOfSample, const std::vector<int> & iMaterials,
                                             std::vector<float> &oStresses, std::vector<float> &oStrains)
{
  oStresses.clear();
  oStrains.clear();

  std::vector<Vector3f> xinit;
  int isample;
  iNumberOfSample = 10;
  for (isample=0; isample<iNumberOfSample; isample++)
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
    computeForceDeformationSample(em, iStepperType, 100);
    float strain = computeStrain(em, iForceAxis);
    xinit = em->x;

    writeVector2File(em->me, "Output//materials.m");
    std::cout << "Strain = " << strain << std::endl;

    /*World * world = new World();
    world->em.push_back(em);

    Render render;
    render.init(world);
    render.loop();*/ 

    delete em;

    oStresses.push_back(forceMagnitude);
    oStrains.push_back(strain);
  }
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

  oMaterials.clear();
  oMaterials.resize(iN[0]*iN[1]*iN[2], indMat2);

  int n[3] = {iN[0], iN[1], iN[2]};
  int axis = idx / (nMat*nMat);
  n[axis] /= 2;

  int i,j,k;
  for (i=0; i<n[0]; i++)
  {
    for (j=0; j<n[1]; j++)
    {
      for (k=0; k<n[2]; k++)
      {
        int voxelIndex = getVoxelIndex(i,j,k, iN[0], iN[1], iN[2]);
        oMaterials[voxelIndex] = indMat1;
      }
    }
  }
}

void SamplesGeneratorImpl::int2CombRandom(std::vector<int> &oMaterials, int idx, int nMat, int nVox)
{
  oMaterials.resize(nVox);
  for(int ii = nVox-1;ii>=0;ii--)
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

void SamplesGeneratorImpl::getMaterialAssignment(int nx, int ny, int nz, const std::vector<int> &iMaterialCombIndices, int nX, int nY, int nZ, const std::vector<std::vector<int> > &iBaseCellMaterials, std::vector<int> &oMaterials)
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

int SamplesGeneratorImpl::computeMaterialParameters(std::vector<MaterialQuad> &iMaterials, std::string iStepperType , 
                                                    const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<std::vector<int> > &iBaseMaterialStructures, int iBaseMatStructureSize[3], int iLevel, int iBlockRep)
{
  int icomb, ncomb=(int)iMaterialAssignments.size();
  for (icomb=0; icomb<ncomb; icomb++)
  {
    int c = 2*iBlockRep;
    int n[3] = {c*iBaseMatStructureSize[0], c*iBaseMatStructureSize[1], c*iBaseMatStructureSize[2]};

    std::vector<int> cellMaterials;
    getMaterialAssignment(n[0], n[1], n[2], iMaterialAssignments[icomb], iBaseMatStructureSize[0], iBaseMatStructureSize[1], iBaseMatStructureSize[2], iBaseMaterialStructures, cellMaterials);

    int iaxis=0;
    for (iaxis=0; iaxis<1; iaxis++)
    {
      //float forceMagnitude = 1.f;
      float forceMagnitude = 0.3f;
      int NumberOfSample = 10;
      int axis = iaxis;
      //int axis = 0;
      Vector3f forceAxis;
      forceAxis[axis] = 1.f;

      std::string strAxis =  (iaxis==0?"x//": "y//");
      std::string FileName = "Output//" + strAxis  + "StressStrain_" + std::to_string(iLevel) + '_' + std::to_string(3) + "_" + std::to_string(icomb) + ".txt";

      std::vector<float> stresses, strains;
      sampleDeformation(n, iMaterials, iStepperType, forceMagnitude, axis, NumberOfSample, cellMaterials, stresses, strains);
      writeData(FileName, forceAxis, iMaterialAssignments[icomb], stresses, strains);
    }
  }
  //exit(0);
  return 0;
}

int SamplesGeneratorImpl::computeMaterialParametersLevel1(std::vector<MaterialQuad> & iMaterials, std::string & iStepperType)
{
  // level 0
  std::vector<std::vector<int> > baseMaterialStructures;
  baseMaterialStructures.resize(2);
  baseMaterialStructures[0].push_back(0);
  baseMaterialStructures[1].push_back(1);
  int N[3] = {1,1,1};

  // level 1
  int level = 1;
  int nVoxCoarse = 8;
  int nmat = 2;
  int ncomb=pow(nmat,nVoxCoarse);
  std::vector<std::vector<int> > materialAssignments(ncomb);
  int icomb;
  for (icomb=0; icomb<ncomb; icomb++)
  {
    int2Comb(materialAssignments[icomb], icomb, nmat, nVoxCoarse);
  }

  std::string FileNameMaterialsBaseStructure = "Output//Materials_" + std::to_string(level)  + ".txt";
  writeMaterialCombinations(FileNameMaterialsBaseStructure,  baseMaterialStructures);

  int blockRep = 3;
  computeMaterialParameters(iMaterials, iStepperType, materialAssignments, baseMaterialStructures, N, level, blockRep);

  return 0;
}

int SamplesGeneratorImpl::computeMaterialParametersLevel2(std::vector<MaterialQuad> & iMaterials, std::string & iStepperType)
{
  std::string fileDirectory = "C://Users//melina//Documents//Projects//fem//vc//runtime//Output//";

  int level = 2;
  int prevLevel = 1;
  int blockSize = 3;

  std::string materialFile = fileDirectory + "Materials_" + std::to_string(prevLevel) + ".txt";
  std::vector<std::vector<int> > baseMaterialStructures;
  readMaterialCombinations(materialFile, baseMaterialStructures);

  std::vector<Scalar> physicalParameters; //YoungModulus, density;
  std::vector<std::vector<int> > materialAssignments;
  int icomb, ncomb=256;
  materialAssignments.resize(ncomb);
  int naxis = 2;
  for (icomb=0; icomb<ncomb; icomb++)
  {
    int iaxis;
    for (iaxis=0; iaxis<naxis; iaxis++)
    {
      std::string subDirectory = (iaxis==0? "x_level1//": "y_level1//");
      std::string sampleFile = fileDirectory + subDirectory + "StressStrain_" + std::to_string(prevLevel) + "_" + std::to_string(blockSize) + "_" + std::to_string(icomb) + ".txt"; 
      Vector3f forceAxis;
      std::vector<int> materials;
      std::vector<Scalar> stresses, strains;
      readData(sampleFile, forceAxis, materials, stresses, strains);

      if (iaxis==0)
      {
        Scalar r=0;
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
          materialAssignments[icomb].insert(materialAssignments[icomb].end(), baseMaterialStructure.begin(), baseMaterialStructure.end());
        }
        int ncell = baseMaterialStructures[0].size()*nmat;
        r /= (Scalar)ncell;
        physicalParameters.push_back(r);
      }
      Scalar Y = computeYoungModulus(strains, stresses);
      physicalParameters.push_back(Y);
    }
  }
  int nparam = naxis + 1;
  std::vector<int> convexHullVertices;
  computeConvexHull(physicalParameters, nparam, convexHullVertices);
  if (1)
  {
    writeVector2File(convexHullVertices, "Output//convexHull.m");
  }
  std::vector<std::vector<int> > newBaseMaterialStructures = getSubVector(materialAssignments, convexHullVertices);

  std::string FileNameMaterials = "Output//Materials_" + std::to_string(level)  + ".txt";
  writeMaterialCombinations(FileNameMaterials,  newBaseMaterialStructures);

  int N[3] = {2,2,2};

  int nVoxCoarse = 8;
  int nmat = (int)convexHullVertices.size();
  //int nNewComb=pow(nmat,nVoxCoarse);
  int nNewComb = 3*nmat*nmat-2;
  std::vector<std::vector<int> > newMaterialAssignments(nNewComb);
  for (icomb=0; icomb<nNewComb; icomb++)
  {
    //int2Comb(newMaterialAssignments[icomb], icomb, nmat, nVoxCoarse);
    int2Comb_biMaterialStructures(newMaterialAssignments[icomb], icomb, nmat, N);
  }
  int blockRep = 3;
  computeMaterialParameters(iMaterials, iStepperType, newMaterialAssignments, newBaseMaterialStructures, N, level, blockRep);

  return 0;
}

int SamplesGeneratorImpl::computeMaterialParameters(std::vector<MaterialQuad> & iMaterials, std::string & iStepperType, int iLevel, int iNbCombPrevLevel)
{
  assert(iLevel>=2);

  std::string fileDirectory = "C://Users//melina//Documents//Projects//fem//vc//runtime//Output//";

  int prevLevel = iLevel-1;
  int blockSize = 3;

  std::string materialFile = fileDirectory + "Materials_" + std::to_string(prevLevel) + ".txt";
  std::vector<std::vector<int> > baseMaterialStructures;
  readMaterialCombinations(materialFile, baseMaterialStructures);

  std::string subDirectories[2];
  subDirectories[0] =  "x_level" + std::to_string(prevLevel) + "//";
  subDirectories[1] =  "y_level" + std::to_string(prevLevel) + "//";

  std::vector<Scalar> physicalParameters; //YoungModulus, density;
  std::vector<std::vector<int> > materialAssignments;
  int icomb, ncomb=iNbCombPrevLevel; //256;
  materialAssignments.resize(ncomb);
  int naxis = 2;
  for (icomb=0; icomb<ncomb; icomb++)
  {
    int iaxis;
    for (iaxis=0; iaxis<naxis; iaxis++)
    {
      std::string sampleFile = fileDirectory + subDirectories[iaxis] + "StressStrain_" + std::to_string(prevLevel) + "_" + std::to_string(blockSize) + "_" + std::to_string(icomb) + ".txt"; 
      Vector3f forceAxis;
      std::vector<int> materials;
      std::vector<Scalar> stresses, strains;
      readData(sampleFile, forceAxis, materials, stresses, strains);

      if (iaxis==0)
      {
        Scalar r=0;
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
          materialAssignments[icomb].insert(materialAssignments[icomb].end(), baseMaterialStructure.begin(), baseMaterialStructure.end());
        }
        int ncell = baseMaterialStructures[0].size()*nmat;
        r /= (Scalar)ncell;
        physicalParameters.push_back(r);
      }
      Scalar Y = computeYoungModulus(strains, stresses);
      physicalParameters.push_back(Y);
    }
  }
  int nparam = naxis + 1;
  std::vector<int> convexHullVertices;
  computeConvexHull(physicalParameters, nparam, convexHullVertices);
  if (1)
  {
    writeVector2File(convexHullVertices, "Output//convexHull.m");
  }
  std::vector<std::vector<int> > newBaseMaterialStructures = getSubVector(materialAssignments, convexHullVertices);

  std::string FileNameMaterials = "Output//Materials_" + std::to_string(iLevel)  + ".txt";
  writeMaterialCombinations(FileNameMaterials,  newBaseMaterialStructures);

  int sideSize = pow(2, prevLevel);
  int N[3] = {sideSize,sideSize,sideSize};

  int nVoxCoarse = 8;
  int nmat = (int)convexHullVertices.size();
  //int nNewComb=pow(nmat,nVoxCoarse);
  int nNewComb = 3*nmat*nmat-2;
  std::vector<std::vector<int> > newMaterialAssignments(nNewComb);
  for (icomb=0; icomb<nNewComb; icomb++)
  {
    //int2Comb(newMaterialAssignments[icomb], icomb, nmat, nVoxCoarse);
    int2Comb_biMaterialStructures(newMaterialAssignments[icomb], icomb, nmat, N);
  }
  int blockRep = 3;
  computeMaterialParameters(iMaterials, iStepperType, newMaterialAssignments, newBaseMaterialStructures, N, iLevel, blockRep);

  return 0;
}



int SamplesGeneratorImpl::run()
{
  std::vector<StrainEneNeo> ene(2);
  ene[0].param[0] = 1;
  ene[0].param[1] = 10;
  ene[1].param[0] = 10;
  ene[1].param[1] = 100;
  std::vector<MaterialQuad> material(ene.size());
  for (unsigned int ii = 0; ii < material.size(); ii++){
    for (unsigned int jj = 0; jj < material[ii].e.size(); jj++){
      material[ii].e[jj] = &ene[ii];
    }
  }
  //std::string stepperType = "newtonCuda";
  std::string stepperType = "newton";

  computeMaterialParametersLevel1(material, stepperType);
  //computeMaterialParametersLevel2(material, stepperType);
  //computeMaterialParameters(material, stepperType, 3, 429);

  return 0 ;

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

  srand(0);

  int shiftRand = 600;
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

  int iHLevel;
  for (iHLevel=MinHierarchyLevel; iHLevel<=MaxHierarchyLevel; iHLevel++)
  {
    std::string FileNameMaterials = "Output//Materials_" + std::to_string(iHLevel)  + ".txt";
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

        std::string FileName = "Output//StressStrain_" + std::to_string(iHLevel) + '_' + std::to_string(iReflevel) + "_" + std::to_string(icomb+shiftRand) + ".txt";

        std::vector<float> stresses, strains;
        sampleDeformation(n, material, stepperType, forceMagnitude, axis, NumberOfSample, gridMaterials, stresses, strains);
        writeData(FileName, forceAxis, coarseCellMaterials, stresses, strains);
      }
    }
    materialCombinations = newmaterialCombinations;
  }

  /*World * world = new World();
  world->em.push_back(em);

  Render render;
  render.init(world);
  render.loop();*/ 
  return 0;
}




