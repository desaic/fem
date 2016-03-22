#include "StrengthAnalysis.h"

#include "EigenUtil.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "ElementRegGrid2D.h"
#include "Element2D.h"

#include "MaterialQuad2D.h"
#include "MaterialQuad.hpp"
#include "StrainLin2D.h"
#include "StrainLin.hpp"

#include "ElementHex2D.h"
#include "ElementHex.hpp"

#include "pardiso_sym.hpp"

#include "MeshUtilities.h"
using namespace meshUtil;

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

#include "cfgUtilities.h"
using namespace cfgUtil;

typedef Eigen::Triplet<cfgScalar> TripletS;

StrengthAnalysis::StrengthAnalysis()
{
  m_physicalSystem = NULL;
  m_pardisoState = NULL;
  m_physicalSystem2D = NULL;

  m_nbBlockRep = 1;
  m_nbSubdiv = 1;

  m_dim = 2;

  m_init = false;
}

StrengthAnalysis::~StrengthAnalysis()
{
  SAFE_DELETE(m_physicalSystem);
  SAFE_DELETE(m_physicalSystem2D);
  if (m_pardisoState)
  {
    int nrows = (int)m_I.size() - 1;
    pardisoFree(m_I.data(), m_J.data(), nrows, m_pardisoState);
    delete m_pardisoState; 
  }
}

void StrengthAnalysis::initPardiso(ElementRegGrid2D *iPhysicalSystem)
{
  bool triangular = true;
  m_I.clear();
  m_J.clear();

  assert(iPhysicalSystem);
  m_K0[0] = iPhysicalSystem->m[0]->getStiffness(iPhysicalSystem->e[0], iPhysicalSystem);
  m_K0[1] = iPhysicalSystem->m[1]->getStiffness(iPhysicalSystem->e[0], iPhysicalSystem);
  if (0)
  {
    std::cout << "K0[0] = " << std::endl;
    std::cout << m_K0[0] << std::endl << std::endl;
    std::cout << "K0[1] = " << std::endl;
    std::cout << m_K0[1] << std::endl << std::endl;
  }

  iPhysicalSystem->stiffnessPattern(m_I, m_J, triangular, true, true, true);

  //pardiso uses 1-based
  for (unsigned int ii = 0; ii < m_I.size(); ii++){
    m_I[ii] ++;
  }
  for (unsigned int ii = 0; ii < m_J.size(); ii++){
    m_J[ii] ++;
  }
  SAFE_DELETE(m_pardisoState);
  m_pardisoState = new PardisoState();
  pardisoInit(m_pardisoState);
  pardisoSymbolicFactorize(m_I.data(), m_J.data(), (int)m_I.size()-1, m_pardisoState);
}

void StrengthAnalysis::initPardiso(ElementRegGrid *iPhysicalSystem)
{
  bool triangular = true;
  m_I.clear();
  m_J.clear();

  assert(iPhysicalSystem);
  m_K0[0] = iPhysicalSystem->m[0]->getStiffness(iPhysicalSystem->e[0], iPhysicalSystem);
  m_K0[1] = iPhysicalSystem->m[1]->getStiffness(iPhysicalSystem->e[0], iPhysicalSystem);
  if (0)
  {
    std::cout << "K0[0] = " << std::endl;
    std::cout << m_K0[0] << std::endl << std::endl;
    std::cout << "K0[1] = " << std::endl;
    std::cout << m_K0[1] << std::endl << std::endl;
  }

  iPhysicalSystem->stiffnessPattern(m_I, m_J, triangular, true, true, true);

  //pardiso uses 1-based
  for (unsigned int ii = 0; ii < m_I.size(); ii++){
    m_I[ii] ++;
  }
  for (unsigned int ii = 0; ii < m_J.size(); ii++){
    m_J[ii] ++;
  }
  SAFE_DELETE(m_pardisoState);
  m_pardisoState = new PardisoState();
  pardisoInit(m_pardisoState);
  pardisoSymbolicFactorize(m_I.data(), m_J.data(), (int)m_I.size()-1, m_pardisoState);
}

ElementRegGrid2D * StrengthAnalysis::createPhysicalSystem(int n[2], std::vector<MaterialQuad2D> &iMaterials, const std::vector<int> &iMatAssignments)
{
  ElementRegGrid2D * physicalSystem = new ElementRegGrid2D(n[0],n[1]);

  // Materials
  int imat=0, nmat=(int)iMaterials.size();
  for (imat=0; imat<nmat; imat++)
  {
    physicalSystem->addMaterial(&iMaterials[imat]);
  }
  assert(physicalSystem->me.size()==iMatAssignments.size());
  physicalSystem->me = iMatAssignments;

  return physicalSystem;
}

ElementRegGrid * StrengthAnalysis::createPhysicalSystem(int n[3], std::vector<MaterialQuad> &iMaterials, const std::vector<int> &iMatAssignments)
{
  ElementRegGrid * physicalSystem = new ElementRegGrid(n[0],n[1],n[2]);

  // Materials
  int imat=0, nmat=(int)iMaterials.size();
  for (imat=0; imat<nmat; imat++)
  {
    physicalSystem->addMaterial(&iMaterials[imat]);
  }
  assert(physicalSystem->me.size()==iMatAssignments.size());
  physicalSystem->me = iMatAssignments;

  return physicalSystem;
}

bool StrengthAnalysis::run2D(int N[2], StructureType iType, std::vector<MaterialQuad2D> &iBaseMaterials, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<cfgScalar> &iParameters, const std::vector<cfgScalar> &iTensors, 
                            int iNbBlockRep, int iNbSubdiv, std::vector<cfgScalar> &oStrengths)
{ 
  oStrengths.clear();

  bool resOk = true;

  std::vector<cfgScalar> E_base; 
  int ibase, nbase=(int)iBaseMaterials.size();
  for (ibase=0; ibase<nbase; ibase++)
  {
    cfgScalar E,  nu;
    cfgScalar lambda = ((StrainLin2D*)(iBaseMaterials[ibase].e[0]))->param[1];
    cfgScalar mu = ((StrainLin2D*)(iBaseMaterials[ibase].e[0]))->param[0];
    fromLamesParametersToYoungModulusPoissonRatio(lambda, mu, E, nu);
    E_base.push_back(E);
  }

  int imat, nmat=(int)iMaterialAssignments.size();
  std::cout << "nmat = " << nmat << std::endl;
  for (imat=0; imat<nmat; imat++)
  {
    if (imat%1000==0)
      std::cout << "imat = " << imat << std::endl;
    std::vector<int> cellMaterials;
    repMaterialAssignment(N[0], N[1], iMaterialAssignments[imat], iNbBlockRep, iNbBlockRep, iNbSubdiv, cellMaterials); 

    int n[2] = {N[0]*iNbBlockRep*iNbSubdiv, N[1]*iNbBlockRep*iNbSubdiv};
    ElementRegGrid2D * physicalSystem = createPhysicalSystem(n, iBaseMaterials, cellMaterials);

    if (!m_init)
    {
      initPardiso(physicalSystem);
      m_init = true;
    }

    int indTensor = 6*imat;
    int ind=0;
    MatrixXS C = MatrixXS::Zero(3, 3);
    std::vector<float> Cvalues;
    for (int i=0; i<3; i++)
    {
      for (int j=0; j<3; j++)
      {
        if (j<=i)
        {
          C(i,j) = iTensors[indTensor + ind++];
          C(j,i) = C(i,j);
        }
      }
    }

    std::vector<cfgScalar> YoungsModuliCoarse;
    int indParam=0;
    if (iType==Cubic)
    {
      indParam = 4*imat;
      cfgScalar ECoarse = iParameters[indParam+1];
      YoungsModuliCoarse.push_back(ECoarse);
      YoungsModuliCoarse.push_back(ECoarse);
    }
    else if (iType==Orthotropic)
    {
      indParam = 5*imat;
      cfgScalar ECoarse1 = iParameters[indParam+1];
      YoungsModuliCoarse.push_back(ECoarse1);
      cfgScalar ECoarse2 = iParameters[indParam+2];
      YoungsModuliCoarse.push_back(ECoarse2);
    }
    cfgScalar strength = computeStrength(physicalSystem, C, YoungsModuliCoarse, E_base);
    oStrengths.push_back(strength);

    delete physicalSystem;
  }
  return resOk;
}

bool StrengthAnalysis::run3D(int N[3], StructureType iType, std::vector<MaterialQuad> &iBaseMaterials, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<cfgScalar> &iParameters, const std::vector<cfgScalar> &iTensors, 
                            int iNbBlockRep, int iNbSubdiv, std::vector<cfgScalar> &oStrengths)
{ 
  oStrengths.clear();

  bool resOk = true;

  std::vector<cfgScalar> E_base; 
  int ibase, nbase=(int)iBaseMaterials.size();
  for (ibase=0; ibase<nbase; ibase++)
  {
    cfgScalar E,  nu;
    cfgScalar lambda = ((StrainLin2D*)(iBaseMaterials[ibase].e[0]))->param[1];
    cfgScalar mu = ((StrainLin2D*)(iBaseMaterials[ibase].e[0]))->param[0];
    fromLamesParametersToYoungModulusPoissonRatio(lambda, mu, E, nu);
    E_base.push_back(E);
  }

  int imat, nmat=(int)iMaterialAssignments.size();
  std::cout << "nmat = " << nmat << std::endl;
  for (imat=0; imat<nmat; imat++)
  {
    if (imat%100==0)
      std::cout << "imat = " << imat << std::endl;
    std::vector<int> cellMaterials;
    repMaterialAssignment(N[0], N[1], N[2], iMaterialAssignments[imat], iNbBlockRep, iNbBlockRep, iNbBlockRep, iNbSubdiv, cellMaterials); 

    int n[3] = {N[0]*iNbBlockRep*iNbSubdiv, N[1]*iNbBlockRep*iNbSubdiv, N[2]*iNbBlockRep*iNbSubdiv};
    ElementRegGrid * physicalSystem = createPhysicalSystem(n, iBaseMaterials, cellMaterials);

    if (!m_init)
    {
      initPardiso(physicalSystem);
      m_init = true;
    }

    int indTensor = 21*imat;
    int ind=0;
    MatrixXS C = MatrixXS::Zero(6, 6);
    std::vector<float> Cvalues;
    for (int i=0; i<6; i++)
    {
      for (int j=0; j<6; j++)
      {
        if (j<=i)
        {
          C(i,j) = iTensors[indTensor + ind++];
          C(j,i) = C(i,j);
        }
      }
    }

    std::vector<cfgScalar> YoungsModuliCoarse;
    int indParam=0;
    if (iType==Cubic)
    {
      indParam = 4*imat;
      cfgScalar ECoarse = iParameters[indParam+1];
      YoungsModuliCoarse.push_back(ECoarse);
      YoungsModuliCoarse.push_back(ECoarse);
      YoungsModuliCoarse.push_back(ECoarse);
    }
    else if (iType==Orthotropic)
    {
      assert(0); // Not implemented
    }
    cfgScalar strength = computeStrength(physicalSystem, C, YoungsModuliCoarse, E_base);
    oStrengths.push_back(strength);

    delete physicalSystem;
  }
  return resOk;
}

cfgScalar StrengthAnalysis::computeStrength(ElementRegGrid2D * iPhysicalSystem, MatrixXS &iElasticityTensorCoarseMaterial, const std::vector<cfgScalar>  &iYoungsModuliCoarse, const std::vector<cfgScalar> &iYoungsModulusBaseMaterials)
{
  cfgScalar strength = FLT_MAX;

  int dim = 2;
  
  Vector2S p(0, 0);
  MatrixXS B = iPhysicalSystem->e[0]->getMatrixB(p, iPhysicalSystem->X);

  bool triangular = true;
  std::vector<double> val;
  MatrixXS K;
  getStiffnessSparse(iPhysicalSystem, val, K);

  int nforce = dim*(int)iPhysicalSystem->fe.size();
  int nrows = (int)m_I.size() - 1;
  pardisoNumericalFactorize(m_I.data(), m_J.data(), val.data(), nrows, m_pardisoState);

  cfgScalar forceMagnitude = 1;
  cfgScalar forces[2][2] = { {forceMagnitude,0}, {0,forceMagnitude} };

  int n[2];
  n[0] = iPhysicalSystem->nx;
  n[1] = iPhysicalSystem->ny;

  std::vector<std::vector<double> > externalForces(2);
  for (int iaxis=0; iaxis<dim; iaxis++)
  {
    externalForces[iaxis].resize(nforce);

    cfgScalar thickness = 1.;
    cfgScalar coeff = forceMagnitude * ((cfgScalar)1/(n[(iaxis+1)%2]*thickness))/2;
    Vector2S ff(0,0);
    ff[iaxis] = coeff;

    std::vector<int> sides;
    std::vector<int> signs;
    sides.push_back(0);
    sides.push_back(1);
    signs.push_back(-1);
    signs.push_back(1);

    int iside=0, nside=(int)signs.size();
    for (iside=0; iside<nside; iside++)
    {
      int indSide = sides[iside];
      int sign = signs[iside];
      std::vector<int> elemIndices;
      std::vector<std::vector<int> > fvIndices;
      getSideVertices(2*iaxis+indSide, iPhysicalSystem, elemIndices, fvIndices);

      int ielem=0, nelem=(int)elemIndices.size();
      for (ielem=0; ielem<nelem; ielem++)
      {
        int indElement = elemIndices[ielem];
        int ivertex=0, nvertex=(int)fvIndices[ielem].size();
        for (ivertex=0; ivertex<nvertex; ivertex++)
        {
          int FvIndex = fvIndices[ielem][ivertex];
          int VertexIndex =iPhysicalSystem->e[indElement]->at(FvIndex);
          externalForces[iaxis][2*VertexIndex] += sign*ff[0];
          externalForces[iaxis][2*VertexIndex+1] += sign*ff[1];
        }
      }
    }
  } 

  std::vector<std::vector<cfgScalar> > u;
  for (unsigned int ii = 0; ii < externalForces.size(); ii++)
  {
    externalForces[ii].resize(nrows);
    std::vector<double> disp(nrows, 0.);
    pardisoBackSubstitute(m_I.data(), m_J.data(), val.data(), nrows, disp.data(), externalForces[ii].data(), m_pardisoState);
    disp.resize(nforce);
    u.push_back(convertVec<double, cfgScalar>(disp));
  }

  int nx = iPhysicalSystem->nx;
  int ny = iPhysicalSystem->ny;
  std::vector<int> corners;
  corners.push_back(iPhysicalSystem->GetVertInd(0, 0));
  corners.push_back(iPhysicalSystem->GetVertInd(0, ny));
  corners.push_back(iPhysicalSystem->GetVertInd(nx, 0));
  corners.push_back(iPhysicalSystem->GetVertInd(nx, ny));
  ElementHex2D coarseElem(corners);
  MatrixXS BCoarse = coarseElem.getMatrixB(p, iPhysicalSystem->X);

  cfgScalar limitStrains[2] = {-1, 1};
  cfgScalar coeff = 0;
  cfgScalar Ebase = 1;
  for (int iaxis=0; iaxis<dim; iaxis++)
  {
    cfgScalar maxStrains[2] = {-FLT_MAX, -FLT_MAX};
    int ielem, nelem=(int)iPhysicalSystem->e.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      Element2D * element = iPhysicalSystem->e[ielem];
      int indMat = iPhysicalSystem->me[ielem];
      std::vector<cfgScalar> uElem = getSubVector<cfgScalar>(u[iaxis], dim, element->getNodeIndices());
      MatrixXS strain = B * toMatrixXS(uElem);
      cfgScalar principleStrains[2];
      getEigenValues2x2SymmetricMatrix(strain.data(), principleStrains);

      cfgScalar maxStrain = principleStrains[0];
      if (maxStrain > maxStrains[indMat])
      {
        maxStrains[indMat] = maxStrain;
      }
    }
    if (maxStrains[0] > 1.e-6)
    {
      coeff = (limitStrains[0]>0? limitStrains[0]/maxStrains[0] : -1);
      Ebase = iYoungsModulusBaseMaterials[0];
    }
    if (maxStrains[1] > 1.e-6)
    {
      if (coeff<=0 || (limitStrains[1]>0 && limitStrains[1]/maxStrains[1] < coeff) )
      {
        coeff = (limitStrains[1]>0? limitStrains[1]/maxStrains[1] : -1);
        Ebase = iYoungsModulusBaseMaterials[1];
      }
    }
    if (coeff > 0)
    {
      std::vector<cfgScalar> uCoarse = getSubVector<cfgScalar>(u[iaxis], dim, corners);
      MatrixXS matUCoarse = toMatrixXS(uCoarse);
      MatrixXS strainCoarse = coeff*BCoarse*matUCoarse;
      MatrixXS stressCoarse = iElasticityTensorCoarseMaterial*strainCoarse/Ebase;

      MatrixXS stressCoarse2 = iYoungsModuliCoarse[iaxis]*strainCoarse/Ebase;
      cfgScalar strength2 = stressCoarse2(iaxis,0);

      std::vector<cfgScalar> stresses;
      stresses.push_back(stressCoarse(0,0));
      stresses.push_back(stressCoarse(1,0));
      stresses.push_back(stressCoarse(2,0));
      strength = std::min(strength, strength2);
      //strength = std::min(strength, cfgMaterialUtilities::computeVonMisesStress2D(stresses));
    }
  }
  return strength;
}

cfgScalar StrengthAnalysis::computeStrength(ElementRegGrid* iPhysicalSystem, MatrixXS &iElasticityTensorCoarseMaterial, const std::vector<cfgScalar>  &iYoungsModuliCoarse, const std::vector<cfgScalar> &iYoungsModulusBaseMaterials)
{
  cfgScalar strength = FLT_MAX;

  int dim = 3;
  
  Vector3S p(0, 0, 0);
  MatrixXS B = ((ElementHex*)(iPhysicalSystem->e[0]))->BMatrix(p, iPhysicalSystem->X);

  bool triangular = true;
  std::vector<double> val;
  MatrixXS K;
  getStiffnessSparse(iPhysicalSystem, val, K);

  int nforce = dim*(int)iPhysicalSystem->fe.size();
  int nrows = (int)m_I.size() - 1;
  pardisoNumericalFactorize(m_I.data(), m_J.data(), val.data(), nrows, m_pardisoState);

  cfgScalar forceMagnitude = 1;
  cfgScalar forces[3][3] = { {forceMagnitude,0,0}, {0,forceMagnitude,0}, {0,0,forceMagnitude}};

  int n[3];
  n[0] = iPhysicalSystem->nx;
  n[1] = iPhysicalSystem->ny;
  n[2] = iPhysicalSystem->nz;

  std::vector<std::vector<double> > externalForces(3);
  for (int iaxis=0; iaxis<dim; iaxis++)
  {
    externalForces[iaxis].resize(nforce);

    cfgScalar coeff = forceMagnitude * ((cfgScalar)1/(n[(iaxis+1)%3]*n[(iaxis+2)%3]))/4;
    Vector3S ff(0,0,0);
    ff[iaxis] = coeff;

    std::vector<int> sides;
    std::vector<int> signs;
    sides.push_back(0);
    sides.push_back(1);
    signs.push_back(-1);
    signs.push_back(1);

    int iside=0, nside=(int)signs.size();
    for (iside=0; iside<nside; iside++)
    {
      int indSide = sides[iside];
      int sign = signs[iside];
      std::vector<int> elemIndices;
      std::vector<std::vector<int> > fvIndices;
      getSideVertices(2*iaxis+indSide, iPhysicalSystem, elemIndices, fvIndices);

      int ielem=0, nelem=(int)elemIndices.size();
      for (ielem=0; ielem<nelem; ielem++)
      {
        int indElement = elemIndices[ielem];
        int ivertex=0, nvertex=(int)fvIndices[ielem].size();
        for (ivertex=0; ivertex<nvertex; ivertex++)
        {
          int FvIndex = fvIndices[ielem][ivertex];
          int VertexIndex =iPhysicalSystem->e[indElement]->at(FvIndex);
          externalForces[iaxis][3*VertexIndex] += sign*ff[0];
          externalForces[iaxis][3*VertexIndex+1] += sign*ff[1];
          externalForces[iaxis][3*VertexIndex+2] += sign*ff[2];
        }
      }
    }
  } 

  std::vector<std::vector<cfgScalar> > u;
  for (unsigned int ii = 0; ii < externalForces.size(); ii++)
  {
    externalForces[ii].resize(nrows);
    std::vector<double> disp(nrows, 0.);
    pardisoBackSubstitute(m_I.data(), m_J.data(), val.data(), nrows, disp.data(), externalForces[ii].data(), m_pardisoState);
    disp.resize(nforce);
    u.push_back(convertVec<double, cfgScalar>(disp));
  }

  int nx = iPhysicalSystem->nx;
  int ny = iPhysicalSystem->ny;
  int nz = iPhysicalSystem->nz;
  std::vector<int> corners;
  corners.push_back(iPhysicalSystem->GetVertInd(0, 0, 0));
  corners.push_back(iPhysicalSystem->GetVertInd(0, 0, nz));
  corners.push_back(iPhysicalSystem->GetVertInd(0, ny, 0));
  corners.push_back(iPhysicalSystem->GetVertInd(0, ny, nz));
  corners.push_back(iPhysicalSystem->GetVertInd(nx, 0, 0));
  corners.push_back(iPhysicalSystem->GetVertInd(nx, 0, nz));
  corners.push_back(iPhysicalSystem->GetVertInd(nx, ny, 0));
  corners.push_back(iPhysicalSystem->GetVertInd(nx, ny, nz));
  ElementHex coarseElem(corners);
  MatrixXS BCoarse = coarseElem.BMatrix(p, iPhysicalSystem->X);

  cfgScalar limitStrains[2] = {-1, 1};
  cfgScalar coeff = 0;
  cfgScalar Ebase = 1;
  for (int iaxis=0; iaxis<dim; iaxis++)
  {
    cfgScalar maxStrains[2] = {-FLT_MAX, -FLT_MAX};
    int ielem, nelem=(int)iPhysicalSystem->e.size();
    for (ielem=0; ielem<nelem; ielem++)
    {
      Element * element = iPhysicalSystem->e[ielem];
      int indMat = iPhysicalSystem->me[ielem];
      std::vector<cfgScalar> uElem = getSubVector<cfgScalar>(u[iaxis], dim, element->getNodeIndices());
      MatrixXS strain = B * toMatrixXS(uElem);
      cfgScalar principleStrains[3];
      getEigenValues3x3SymmetricMatrix(strain.data(), principleStrains);

      cfgScalar maxStrain = principleStrains[0];
      if (maxStrain > maxStrains[indMat])
      {
        maxStrains[indMat] = maxStrain;
      }
    }
    if (maxStrains[0] > 1.e-6)
    {
      coeff = (limitStrains[0]>0? limitStrains[0]/maxStrains[0] : -1);
      Ebase = iYoungsModulusBaseMaterials[0];
    }
    if (maxStrains[1] > 1.e-6)
    {
      if (coeff<=0 || (limitStrains[1]>0 && limitStrains[1]/maxStrains[1] < coeff) )
      {
        coeff = (limitStrains[1]>0? limitStrains[1]/maxStrains[1] : -1);
        Ebase = iYoungsModulusBaseMaterials[1];
      }
    }
    if (coeff > 0)
    {
      std::vector<cfgScalar> uCoarse = getSubVector<cfgScalar>(u[iaxis], dim, corners);
      MatrixXS matUCoarse = toMatrixXS(uCoarse);
      MatrixXS strainCoarse = coeff*BCoarse*matUCoarse;
      MatrixXS stressCoarse = iElasticityTensorCoarseMaterial*strainCoarse/Ebase;

      MatrixXS stressCoarse2 = iYoungsModuliCoarse[iaxis]*strainCoarse/Ebase;
      cfgScalar strength2 = stressCoarse2(iaxis,0);

      std::vector<cfgScalar> stresses;
      stresses.push_back(stressCoarse(0,0));
      stresses.push_back(stressCoarse(1,0));
      stresses.push_back(stressCoarse(2,0));
      strength = std::min(strength, strength2);
      //strength = std::min(strength, cfgMaterialUtilities::computeVonMisesStress2D(stresses));
    }
  }
  return strength;
}

void StrengthAnalysis::getStiffnessSparse(ElementRegGrid2D * iPhysicalSystem, std::vector<double> &oValues, MatrixXS &oK)
{
  bool fixRigid = true;
  bool periodic = true;
  bool constrained = false;
  bool trig = true;
  int dim = 2;
  int N = dim * (int)iPhysicalSystem->x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N, N);
  for (unsigned int ii = 0; ii<iPhysicalSystem->e.size(); ii++){
    Element2D * ele = iPhysicalSystem->e[ii];
    int nV = ele->nV();
    MatrixXS &K = m_K0[iPhysicalSystem->me[ii]];

    for (int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for (int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for (int dim1 = 0; dim1<dim; dim1++){
          for (int dim2 = 0; dim2<dim; dim2++){
            if (trig && (dim * vk + dim2 > dim * vj + dim1)) {
              continue;
            }
            cfgScalar val = K(dim * jj + dim1, dim * kk + dim2);
            if (constrained){
              if ( (iPhysicalSystem->fixed[vk] || iPhysicalSystem->fixed[vj]) 
                && (vj != vk || dim1 != dim2)){
                val = 0;
              }
            }
            if (dim * vj + dim1 == dim * vk + dim2){
              val *= 1 + (cfgScalar)1e-5;
            }
            TripletS triple(dim * vj + dim1, dim * vk + dim2, val);
            coef.push_back(triple);
          }
        }
      }
    }
  }
  Ksparse.setFromTriplets(coef.begin(), coef.end());
  oK = Ksparse;
  if (fixRigid)
  {
    iPhysicalSystem->fixTranslation(Ksparse, trig, iPhysicalSystem);
    iPhysicalSystem->fixRotation(Ksparse, trig, iPhysicalSystem);
  }
  if (periodic)
  {
    iPhysicalSystem->enforcePeriodicity(Ksparse, trig, iPhysicalSystem);
  }

  oValues.resize(Ksparse.nonZeros());
  int idx = 0;
  for (int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<cfgScalar>::InnerIterator it(Ksparse, ii); it; ++it){
      oValues[idx] = it.value();
      idx++;
    }
  }
}

void StrengthAnalysis::getStiffnessSparse(ElementRegGrid* iPhysicalSystem, std::vector<double> &oValues, MatrixXS &oK)
{
  bool fixRigid = true;
  bool periodic = true;
  bool constrained = false;
  bool trig = true;
  int dim = 3;
  int N = dim * (int)iPhysicalSystem->x.size();
  std::vector<TripletS> coef;
  Eigen::SparseMatrix<cfgScalar> Ksparse(N, N);
  for (unsigned int ii = 0; ii<iPhysicalSystem->e.size(); ii++){
    Element * ele = iPhysicalSystem->e[ii];
    int nV = ele->nV();
    MatrixXS &K = m_K0[iPhysicalSystem->me[ii]];

    for (int jj = 0; jj<nV; jj++){
      int vj = ele->at(jj);
      for (int kk = 0; kk<nV; kk++){
        int vk = ele->at(kk);
        for (int dim1 = 0; dim1<dim; dim1++){
          for (int dim2 = 0; dim2<dim; dim2++){
            if (trig && (dim * vk + dim2 > dim * vj + dim1)) {
              continue;
            }
            cfgScalar val = K(dim * jj + dim1, dim * kk + dim2);
            if (constrained){
              if ( (iPhysicalSystem->fixed[vk] || iPhysicalSystem->fixed[vj]) 
                && (vj != vk || dim1 != dim2)){
                val = 0;
              }
            }
            if (dim * vj + dim1 == dim * vk + dim2){
              val *= 1 + (cfgScalar)1e-5;
            }
            TripletS triple(dim * vj + dim1, dim * vk + dim2, val);
            coef.push_back(triple);
          }
        }
      }
    }
  }
  Ksparse.setFromTriplets(coef.begin(), coef.end());
  oK = Ksparse;
  if (fixRigid)
  {
    iPhysicalSystem->fixTranslation(Ksparse, trig, iPhysicalSystem);
    iPhysicalSystem->fixRotation(Ksparse, trig, iPhysicalSystem);
  }
  if (periodic)
  {
    iPhysicalSystem->enforcePeriodicity(Ksparse, trig, iPhysicalSystem);
  }

  oValues.resize(Ksparse.nonZeros());
  int idx = 0;
  for (int ii = 0; ii<Ksparse.rows(); ii++){
    for (Eigen::SparseMatrix<cfgScalar>::InnerIterator it(Ksparse, ii); it; ++it){
      oValues[idx] = it.value();
      idx++;
    }
  }
}


