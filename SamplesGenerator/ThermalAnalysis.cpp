#include "ThermalAnalysis.h"

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

ThermalAnalysis::ThermalAnalysis()
{
  m_physicalSystem = NULL;
  m_pardisoState = NULL;
  m_physicalSystem2D = NULL;

  m_nbBlockRep = 1;
  m_nbSubdiv = 1;

  m_dim = 2;

  m_init = false;
}

ThermalAnalysis::~ThermalAnalysis()
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

void ThermalAnalysis::initPardiso(ElementRegGrid2D *iPhysicalSystem)
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

void ThermalAnalysis::initPardiso(ElementRegGrid *iPhysicalSystem)
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

ElementRegGrid2D * ThermalAnalysis::createPhysicalSystem(int n[2], std::vector<MaterialQuad2D> &iMaterials, const std::vector<int> &iMatAssignments)
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

ElementRegGrid * ThermalAnalysis::createPhysicalSystem(int n[3], std::vector<MaterialQuad> &iMaterials, const std::vector<int> &iMatAssignments)
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

bool ThermalAnalysis::run2D(int N[2], std::vector<MaterialQuad2D> &iBaseMaterials, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<cfgScalar> &iParameters, const std::vector<cfgScalar> &iTensors, 
                            int iNbBlockRep, int iNbSubdiv, std::vector<cfgScalar> &oThermalCoefficients)
{ 
  oThermalCoefficients.clear();

  bool resOk = true;

  cfgScalar E,  nu;
  cfgScalar deltaT = 1;
  cfgScalar alpha = 1;
  cfgScalar lambda = ((StrainLin2D*)(iBaseMaterials[1].e[0]))->param[1];
  cfgScalar mu = ((StrainLin2D*)(iBaseMaterials[1].e[0]))->param[0];
  fromLamesParametersToYoungModulusPoissonRatio(lambda, mu, E, nu);

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
    int indParam = 4*imat;
    cfgScalar ECoarse = iParameters[indParam+1];
    cfgScalar nuCoarse = iParameters[indParam+2];
    cfgScalar homogenizedThermalCoefficient = computeThermalCoefficient(physicalSystem, E, nu, deltaT, alpha, C, ECoarse, nuCoarse);
    oThermalCoefficients.push_back(homogenizedThermalCoefficient);

    delete physicalSystem;
  }
  return resOk;
}

bool ThermalAnalysis::run3D(int N[3], std::vector<MaterialQuad> &iBaseMaterials, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<cfgScalar> &iParameters, const std::vector<cfgScalar> &iTensors, 
                            int iNbBlockRep, int iNbSubdiv, std::vector<cfgScalar> &oThermalCoefficients)
{ 
  oThermalCoefficients.clear();

  bool resOk = true;

  cfgScalar E,  nu;
  cfgScalar deltaT = 1;
  cfgScalar alpha = 1;
  cfgScalar lambda = ((StrainLin2D*)(iBaseMaterials[1].e[0]))->param[1];
  cfgScalar mu = ((StrainLin2D*)(iBaseMaterials[1].e[0]))->param[0];
  fromLamesParametersToYoungModulusPoissonRatio(lambda, mu, E, nu);

  int imat, nmat=(int)iMaterialAssignments.size();
  std::cout << "nmat = " << nmat << std::endl;
  for (imat=0; imat<nmat; imat++)
  {
    //if (imat%100==0)
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
    int indParam = 4*imat;
    cfgScalar ECoarse = iParameters[indParam+1];
    cfgScalar nuCoarse = iParameters[indParam+2];
    cfgScalar homogenizedThermalCoefficient = computeThermalCoefficient(physicalSystem, E, nu, deltaT, alpha, C, ECoarse, nuCoarse);

    cfgScalar density = iParameters[indParam];
    cfgScalar K1 = E/(3*(1-2*nu));
    cfgScalar K2 = 1./1000 * E/(3*(1-2*0.1));
    cfgScalar K = ECoarse/(3*(1-2*nuCoarse));
    cfgScalar val = density + 1/(1/K1 - 1/K2) * (1/K - (density/K1 + (1-density)/K2));

    //oThermalCoefficients.push_back(val);
    //oThermalCoefficients.push_back(homogenizedThermalCoefficient);

    delete physicalSystem;
  }
  return resOk;
}

cfgScalar ThermalAnalysis::computeThermalCoefficient(ElementRegGrid2D * iPhysicalSystem, cfgScalar iYoungModulus, cfgScalar iPoissonRatio, cfgScalar iDeltaT, cfgScalar iCoeffThermalExpansion, 
                                                     MatrixXS &iElasticityTensor, cfgScalar iYoungModulusCoarse, cfgScalar iPoissonRatioCoarse)
{
  cfgScalar coeff = 0.5*iYoungModulus*iDeltaT/(1-iPoissonRatio);
  coeff = 1;

  int dim = 2;
  std::vector<double> heatForces(dim*iPhysicalSystem->fe.size());

  Vector2S p(0, 0);
  MatrixXS B = iPhysicalSystem->e[0]->getMatrixB(p, iPhysicalSystem->X);
  MatrixXS It = MatrixXS::Zero(1, 3);
  It(0,0) = 1;
  It(0,1) = 1;
  MatrixXS ItB = It*B;
  MatrixXS coeffItB = coeff*iCoeffThermalExpansion*ItB;

  int ielem, nelem=(int)iPhysicalSystem->e.size();
  for (ielem=0; ielem<nelem; ielem++)
  {
    Element2D * element = iPhysicalSystem->e[ielem];
    int indMat = iPhysicalSystem->me[ielem];
    if (indMat == 1)
    {
      int ivertex, nvertex = element->nV();
      for (ivertex=0; ivertex<nvertex; ivertex++)
      {
        int indVertex = element->at(ivertex);
        for (int icoord=0; icoord<dim; icoord++)
        {
          heatForces[dim*indVertex + icoord] += coeffItB(0, dim*ivertex + icoord);
        }
      }
    }
  }

  bool triangular = true;
  std::vector<double> val;
  MatrixXS K;
  getStiffnessSparse(iPhysicalSystem, val, K);
  std::vector<float> Kvalues;

  int nforce = heatForces.size();
  int nrows = (int)m_I.size() - 1;
  pardisoNumericalFactorize(m_I.data(), m_J.data(), val.data(), nrows, m_pardisoState);

  std::vector<double> externalForces = heatForces;
  externalForces.resize(nrows);

  std::vector<double> u(nrows, 0.);
  pardisoBackSubstitute(m_I.data(), m_J.data(), val.data(), nrows, u.data(), externalForces.data(), m_pardisoState);
  u.resize(nforce);

  cfgScalar c = iCoeffThermalExpansion*iYoungModulus/(1-iPoissonRatio);
  cfgScalar Wel = innerProd(heatForces, u);
  cfgScalar Wth = 0;
  for (ielem=0; ielem<nelem; ielem++)
  {
    Element2D * element = iPhysicalSystem->e[ielem];
    int indMat = iPhysicalSystem->me[ielem];
    if (indMat == 1)
    {
      std::vector<cfgScalar> uElem = convertVec<double, cfgScalar>(getSubVector<double>(u, dim, element->getNodeIndices()));
      MatrixXS strain = B * toMatrixXS(uElem);
      Wth += (strain(0,0)+strain(1,0));
    }
  }
  Wth *= c/(float)(iPhysicalSystem->nx * iPhysicalSystem->ny);

  int nx = iPhysicalSystem->nx;
  int ny = iPhysicalSystem->ny;
  std::vector<int> corners;
  corners.push_back(iPhysicalSystem->GetVertInd(0, 0));
  corners.push_back(iPhysicalSystem->GetVertInd(0, ny));
  corners.push_back(iPhysicalSystem->GetVertInd(nx, 0));
  corners.push_back(iPhysicalSystem->GetVertInd(nx, ny));
  ElementHex2D coarseElem(corners);

  MatrixXS BCoarse = coarseElem.getMatrixB(p, iPhysicalSystem->X);
  std::vector<cfgScalar> uCoarse = convertVec<double, cfgScalar>(getSubVector<double>(u, dim, corners));
  MatrixXS matUCoarse = toMatrixXS(uCoarse);

  MatrixXS strainCoarse = BCoarse*matUCoarse;
  MatrixXS ItStrainCoarse = It*strainCoarse;

  cfgScalar homogenizedThermalCoefficient = 0;
  cfgScalar ItStrain = ItStrainCoarse(0,0);
  if (fabs(ItStrain)>1.e-8)
  {
   // homogenizedThermalCoefficient = W/(ItStrain*coeff);
  }

  MatrixXS matU = toMatrixXS(convertVec<double, cfgScalar>(u));
  MatrixXS uKu = matU.transpose()*K*matU;
  //MatrixXS uKu_coarse = matUCoarse.transpose()*KCoarse*matUCoarse;

  MatrixXS KCoarse = BCoarse.transpose()*iElasticityTensor*BCoarse;
  MatrixXS Ku_coarse = KCoarse*matUCoarse;

  cfgScalar E = iYoungModulusCoarse;
  cfgScalar nu = iPoissonRatioCoarse;
  cfgScalar coeffCoarse = 0;
  if (fabs(ItStrain)>1.e-8)
  {
    coeffCoarse = Wth/(strainCoarse(0,0)+strainCoarse(1,0));

    //cfgScalar c = iYoungModulus/(1-2*iPoissonRatio);
    homogenizedThermalCoefficient = coeffCoarse*(1-nu)/E;
  }
  return homogenizedThermalCoefficient;
}

cfgScalar ThermalAnalysis::computeThermalCoefficient(ElementRegGrid * iPhysicalSystem, cfgScalar iYoungModulus, cfgScalar iPoissonRatio, cfgScalar iDeltaT, cfgScalar iCoeffThermalExpansion, 
                                                     MatrixXS &iElasticityTensor, cfgScalar iYoungModulusCoarse, cfgScalar iPoissonRatioCoarse)
{
  cfgScalar coeff = 0.5*iYoungModulus*iDeltaT/(1-2*iPoissonRatio);
  coeff = 1;

  int dim = 3;
  std::vector<double> heatForces(dim*iPhysicalSystem->fe.size());

  Vector3S p(0, 0, 0);
  MatrixXS B = ((ElementHex*)(iPhysicalSystem->e[0]))->BMatrix(p, iPhysicalSystem->X);
  MatrixXS It = MatrixXS::Zero(1, 6);
  It(0,0) = 1;
  It(0,1) = 1;
  It(0,2) = 1;
  MatrixXS ItB = It*B;
  MatrixXS coeffItB = coeff*iCoeffThermalExpansion*ItB;

  int ielem, nelem=(int)iPhysicalSystem->e.size();
  for (ielem=0; ielem<nelem; ielem++)
  {
    Element * element = iPhysicalSystem->e[ielem];
    int indMat = iPhysicalSystem->me[ielem];
    if (indMat == 1)
    {
      int ivertex, nvertex = element->nV();
      for (ivertex=0; ivertex<nvertex; ivertex++)
      {
        int indVertex = element->at(ivertex);
        for (int icoord=0; icoord<dim; icoord++)
        {
          heatForces[dim*indVertex + icoord] += coeffItB(0, dim*ivertex + icoord);
        }
      }
    }
  }

  bool triangular = true;
  std::vector<double> val;
  MatrixXS K;
  getStiffnessSparse(iPhysicalSystem, val, K);
  std::vector<float> Kvalues;

  int nforce = heatForces.size();
  int nrows = (int)m_I.size() - 1;
  pardisoNumericalFactorize(m_I.data(), m_J.data(), val.data(), nrows, m_pardisoState);

  std::vector<double> externalForces = heatForces;
  externalForces.resize(nrows);

  std::vector<double> u(nrows, 0.);
  pardisoBackSubstitute(m_I.data(), m_J.data(), val.data(), nrows, u.data(), externalForces.data(), m_pardisoState);
  u.resize(nforce);

  cfgScalar c = iCoeffThermalExpansion*iYoungModulus/(1-2*iPoissonRatio);
  cfgScalar Wel = innerProd(heatForces, u);
  cfgScalar Wth = 0;
  for (ielem=0; ielem<nelem; ielem++)
  {
    Element * element = iPhysicalSystem->e[ielem];
    int indMat = iPhysicalSystem->me[ielem];
    if (indMat == 1)
    {
      std::vector<cfgScalar> uElem = convertVec<double, cfgScalar>(getSubVector<double>(u, dim, element->getNodeIndices()));
      MatrixXS strain = B * toMatrixXS(uElem);
      Wth += (strain(0,0)+strain(1,0)+strain(2,0));
    }
  }
  Wth *= c/(float)(iPhysicalSystem->nx * iPhysicalSystem->ny * iPhysicalSystem->nz);

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
  std::vector<cfgScalar> uCoarse = convertVec<double, cfgScalar>(getSubVector<double>(u, dim, corners));
  MatrixXS matUCoarse = toMatrixXS(uCoarse);

  MatrixXS strainCoarse = BCoarse*matUCoarse;
  MatrixXS ItStrainCoarse = It*strainCoarse;

  cfgScalar homogenizedThermalCoefficient = 0;
  cfgScalar ItStrain = ItStrainCoarse(0,0);
  if (fabs(ItStrain)>1.e-8)
  {
   // homogenizedThermalCoefficient = W/(ItStrain*coeff);
  }

  MatrixXS matU = toMatrixXS(convertVec<double, cfgScalar>(u));
  MatrixXS uKu = matU.transpose()*K*matU;
  //MatrixXS uKu_coarse = matUCoarse.transpose()*KCoarse*matUCoarse;

  MatrixXS KCoarse = BCoarse.transpose()*iElasticityTensor*BCoarse;
  MatrixXS Ku_coarse = KCoarse*matUCoarse;

  cfgScalar E = iYoungModulusCoarse;
  cfgScalar nu = iPoissonRatioCoarse;
  cfgScalar coeffCoarse = 0;
  if (fabs(ItStrain)>1.e-8)
  {
    coeffCoarse = Wth/(strainCoarse(0,0)+strainCoarse(1,0)+strainCoarse(2,0));

    //cfgScalar c = iYoungModulus/(1-2*iPoissonRatio);
    homogenizedThermalCoefficient = coeffCoarse*(1-2*nu)/E;
  }

  return homogenizedThermalCoefficient;
}

/*
bool ThermalAnalysis::run3D(int N[3], const std::vector<MaterialQuad> &iBaseMaterials, const std::vector<int> &iMaterialAssignments, const std::vector<float> &iCurrentParams, 
                              const std::vector<float> &iTargetParams, std::vector<std::vector<int> > &oNewMaterialAssignments)
{ 
  oNewMaterialAssignments.clear();

  assert(iTargetParams.size()==4);
  float density = iTargetParams[0];
  float E = iTargetParams[1];
  float nu = iTargetParams[2];
  float mu = iTargetParams[3];

  float currentDensity = iCurrentParams[0];
  float currentE = iCurrentParams[1];
  float currentNu = iCurrentParams[2];
  float currentMu = iCurrentParams[3];

  float minDensity = 1.f/(N[0]*N[1]*N[2]);
  float maxDensity = 1-1.f/(N[0]*N[1]*N[2]);
  if (density < minDensity)
  {
    float d = currentDensity-density;
    float ratio = (currentDensity-minDensity)/d;
    density = currentDensity + ratio*(density-currentDensity);
 
    Vector3S p0(currentE, currentNu, currentMu);
    Vector3S p(E, nu, mu);
    if ( (p-p0).squaredNorm() < 1.e-6 && (currentDensity-minDensity)<1./(N[0]*N[1]*N[2]) )
    {
      // too small variation in parameters
      return true; 
    }
    else
    {
      p = p0 + ratio*(p-p0);
      E = p[0];
      nu = p[1];
      mu = p[2];
    }
  }
  else if (density > maxDensity)
  {
    float d = density-currentDensity;
    float ratio = (maxDensity-currentDensity)/d;
    density = currentDensity + ratio*(density-currentDensity);
 
    Vector3S p0(currentE, currentNu, currentMu);
    Vector3S p(E, nu, mu);
    if ( (p-p0).squaredNorm() < 1.e-6 && (maxDensity-currentDensity)<1./(N[0]*N[1]*N[2]) )
    {
      // too small variation in parameters
      return true; 
    }
    else
    {
      p = p0 + ratio*(p-p0);
      E = p[0];
      nu = p[1];
      mu = p[2];
    }
  }
 

  bool resOk = true;

  int nx = N[0], ny = N[1], nz = N[2];

  ElementRegGrid * em = new ElementRegGrid(nx,ny,nz);
  std::vector<StrainLin> ene(1);

  std::vector<MaterialQuad> baseMaterials = iBaseMaterials;
  em->addMaterial(&baseMaterials[1]);
  em->initArrays();
  em->check();

  FEM3DFun * fem = new FEM3DFun();
  RealField * field = NULL;
  if (m_structureType == Cubic)
  {
    PiecewiseConstantCubic3D * cubicfield = new PiecewiseConstantCubic3D();
    cubicfield->allocate(nx/2, ny/2, nz/2);
    field = cubicfield;
  }
  else if (m_structureType == Orthotropic)
  {
    PiecewiseConstantSym3D * symfield = new PiecewiseConstantSym3D();
    symfield->allocate(nx/2, ny/2, nz/2);
    field = symfield;
  }
  else
  {
    PiecewiseConstant3D * constantfield = new PiecewiseConstant3D();
    constantfield->allocate(nx/2, ny/2, nz/2);
    field = constantfield;
  }
  fem->field = field;
  fem->lowerBounds = 1e-3 * Eigen::VectorXd::Ones(field->param.size());
  fem->upperBounds = Eigen::VectorXd::Ones(field->param.size());

  fem->em = em;
  fem->gridSize[0] = nx;
  fem->gridSize[1] = ny;
  fem->gridSize[2] = nz;

  std::vector<int> matAssignment;
  if (m_structureType == Cubic)
  {
    getQuarter(nx, ny, nz ,iMaterialAssignments, matAssignment);
  }
  else
  {
    matAssignment = iMaterialAssignments;
  }

  double shrink = 0.3;
  Eigen::VectorXd x0(field->param.size());
  assert(x0.rows()==(int)matAssignment.size());
  for (int i=0; i<x0.rows(); i++)
  {
    double old_density;
    if (matAssignment[i]==0)
    {
      old_density = 0;
    }
    else
    {
      old_density = 1;
    }
    x0[i] = 0.5 + (old_density - 0.5) * shrink;
  } 
  float fx = 1;
  fem->forceMagnitude = (double)fx;

  fem->init(x0);
  fem->setParam(x0);
  double val = fem->f();

  mu = (float)fem->G(3,3);

  Eigen::MatrixXd target_G = computeTargetStrains(3, E, nu, mu);
  if (1)
  {
    std::cout << "target G = " << std::endl;
    std::cout << target_G << std::endl << std::endl;

    float current_E = iCurrentParams[1];
    float current_nu = iCurrentParams[2];
    float current_mu = (float)fem->G(3,3);
    Eigen::MatrixXd current_G = computeTargetStrains(3, current_E, current_nu, current_mu);

    std::cout << "current G = " << std::endl;
    std::cout << current_G << std::endl << std::endl;
  }

  std::cout << "target density = " << density << std::endl;

  fem->m0 = density; //0.5 * sum(fem->distribution) / fem->distribution.size();
  fem->mw = 0.01*fem->G(0, 0) / fem->density;
  double weightMu = 0.;
  fem->G0 = target_G;
  fem->wG(0) = 1;
  fem->wG(1) = 1;
  fem->wG(2) = 1;
  fem->wG(3) = weightMu;
  fem->wG(4) = weightMu;
  fem->wG(5) = weightMu;

  int nbSteps = 500;
  std::vector<Eigen::VectorXd> parameters;
  gradientDescent(fem, x0, nbSteps, parameters);

  std::set<std::vector<int> > newMaterials;
  int ipoint=0, npoint=(int)parameters.size();
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    const Eigen::VectorXd & x0 = parameters[ipoint];
    std::vector<int> mat(x0.rows());
    for (int i=0; i<x0.rows(); i++)
    {
      if (x0[i] < 0.5)
      {
        mat[i] = 0;
      }
      else
      {
        mat[i] = 1;
      }
    }
    if (newMaterials.count(mat)==0)
    {
      //dumpStructure(nx/2, ny/2, mat);

      newMaterials.insert(mat);
      if (m_structureType == Cubic)
      {
        oNewMaterialAssignments.push_back(mat);
      }
      else
      {
        oNewMaterialAssignments.push_back(mat);
      }
    }
  }
  delete em;
  delete fem;

  return resOk;
}
*/ 

void ThermalAnalysis::getStiffnessSparse(ElementRegGrid2D * iPhysicalSystem, std::vector<double> &oValues, MatrixXS &oK)
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
              val *= 1 + 1e-5;
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

void ThermalAnalysis::getStiffnessSparse(ElementRegGrid* iPhysicalSystem, std::vector<double> &oValues, MatrixXS &oK)
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
              val *= 1 + 1e-5;
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