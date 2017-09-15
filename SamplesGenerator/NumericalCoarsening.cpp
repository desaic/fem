
#include "NumericalCoarsening.h"

#include "ElementRegGrid2D.h"
#include "ElementRegGrid.hpp"
#include "MaterialQuad2D.h"
#include "MaterialQuad.hpp"
#include "Stepper2D.h"
#include "Stepper.hpp"
#include "Element2D.h"
#include "Element.hpp"
#include "StepperNewton2D.h"
#include "StepperNewton.hpp"
#include "pardiso_sym.hpp"
#include "Timer.hpp"
#include "StrainEne2D.h"
#include "StrainEne.hpp"

#include "MeshUtilities.h"
using namespace meshUtil;

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

#include "cfgUtilities.h"
using namespace cfgUtil;

typedef Eigen::Triplet<cfgScalar> TripletS;

NumericalCoarsening::NumericalCoarsening()
{
  m_physicalSystem = NULL;
  m_pardisoState = NULL;
  m_physicalSystem = NULL;
  m_physicalSystem2D = NULL;

  m_nbBlockRep = 1;
  m_nbSubdiv = 1;

  m_dim = 2;

  m_init = false;
}

NumericalCoarsening::~NumericalCoarsening()
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

void NumericalCoarsening::initPardiso()
{
  //Timer t;
  bool triangular = true;
  m_I.clear();
  m_J.clear();

  m_K0.clear();
  if (m_dim==2)
  {
    assert(m_physicalSystem2D);
    int nmat = (int)m_physicalSystem2D->m.size();
    for (int imat=0; imat<nmat; imat++)
    {
      m_K0.push_back(m_physicalSystem2D->m[imat]->getStiffness(m_physicalSystem2D->e[0], m_physicalSystem2D));
    }
    //t.start();
    m_physicalSystem2D->stiffnessPattern(m_I, m_J, triangular, true, true, true);
    //t.end();
    //std::cout << " stiffnessPattern: " <<  t.getSeconds() << std::endl;
  }
  else
  {
    assert(m_physicalSystem);
    int nmat = (int)m_physicalSystem->m.size();
    for (int imat=0; imat<nmat; imat++)
    {
      m_K0.push_back(m_physicalSystem->m[imat]->getStiffness(m_physicalSystem->e[0], m_physicalSystem));
    }
    //t.start();
    m_physicalSystem->stiffnessPattern(m_I, m_J, triangular, true, true, true);
    //t.end();
    //std::cout << " stiffnessPattern: " <<  t.getSeconds() << std::endl;
  }
  if (0)
  {
    int imat=0;
    std::cout << "K0[" << imat << "] = " << std::endl;
    std::cout << m_K0[imat] << std::endl << std::endl;
  }
  //pardiso uses 1-based
  for (unsigned int ii = 0; ii < m_I.size(); ii++){
    m_I[ii] ++;
  }
  for (unsigned int ii = 0; ii < m_J.size(); ii++){
    m_J[ii] ++;
  }
  SAFE_DELETE(m_pardisoState);
  m_pardisoState = new PardisoState();
  //t.start();
  pardisoInit(m_pardisoState);
  //t.end();
  //std::cout << " pardisoInit: " <<  t.getSeconds() << std::endl;
  //t.start();
  pardisoSymbolicFactorize(m_I.data(), m_J.data(), (int)m_I.size()-1, m_pardisoState);
  //t.end();
  //std::cout << " pardisoSymbolicFactorize: " <<  t.getSeconds() << std::endl;
}

void NumericalCoarsening::init(int N[3], int iNbBlockRep, int iNbSubdiv, const std::vector<MaterialQuad> &iBaseMaterials)
{
  m_init = true;

  m_dim = 3;
  m_nbBlockRep = iNbBlockRep;
  m_nbSubdiv = iNbSubdiv;

  int n[3] = {N[0]*m_nbBlockRep*m_nbSubdiv, N[1]*m_nbBlockRep*m_nbSubdiv, N[2]*m_nbBlockRep*m_nbSubdiv};

  m_baseMaterials3D = iBaseMaterials;

  std::vector<int> cellMaterials(n[0]*n[1]*n[2], 0);
  SAFE_DELETE(m_physicalSystem);
  m_physicalSystem = createPhysicalSystem(n, m_baseMaterials3D, cellMaterials);

  initPardiso();
}

void NumericalCoarsening::init(int N[2], int iNbBlockRep, int iNbSubdiv, const std::vector<MaterialQuad2D> &iBaseMaterials)
{
  m_init = true;

  m_dim = 2;
  m_nbBlockRep = iNbBlockRep;
  m_nbSubdiv = iNbSubdiv;

  int n[2] = {N[0]*m_nbBlockRep*m_nbSubdiv, N[1]*m_nbBlockRep*m_nbSubdiv};

  m_baseMaterials2D = iBaseMaterials;

  std::vector<int> cellMaterials(n[0]*n[1], 0);
  SAFE_DELETE(m_physicalSystem);
  m_physicalSystem2D = createPhysicalSystem(n, m_baseMaterials2D, cellMaterials);

  initPardiso();
}

void NumericalCoarsening::computeCoarsenedElasticityTensorAndParameters(const std::vector<int> &iMaterials, int N[2], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad2D> &iBaseMaterials, 
                                                                        StructureType iType, std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters)
{
  std::vector<int> cellMaterials;
  repMaterialAssignment(N[0], N[1], iMaterials, iNbBlockRep, iNbBlockRep, iNbSubdiv, cellMaterials); 

  int n[2] = {N[0]*iNbBlockRep*iNbSubdiv, N[1]*iNbBlockRep*iNbSubdiv};
  ElementRegGrid2D * physicalSystem = createPhysicalSystem(n, iBaseMaterials, cellMaterials);
  Stepper2D * stepper = createStepper(physicalSystem);

  cfgScalar forceMagnitude = 1;
  std::vector<std::vector<cfgScalar> > harmonicDisp;
  computeHarmonicDisplacements(physicalSystem, stepper, forceMagnitude, iType, harmonicDisp);

  std::vector<Vector2S> h = cfgMaterialUtilities::toVector2S(harmonicDisp[0]);
  std::vector<Vector2S> x = cfgUtil::add(physicalSystem->X, h);

  MatrixXS C;
  meshUtil::computeCoarsenedElasticityTensor(*physicalSystem, harmonicDisp, C);

  std::vector<float> Cvalues;
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      if (j<=i)
      {
        ioTensorValues.push_back(C(i,j));
      }
    }
  }

  std::vector<std::vector<cfgScalar> > strains;
  meshUtil::computeStrains(*physicalSystem, harmonicDisp, strains);

  std::vector<std::vector<int> > baseMaterials(2);
  baseMaterials[0].push_back(0);
  baseMaterials[1].push_back(1);

  cfgScalar r = computeMaterialRatio(iMaterials, baseMaterials);
  ioParameters.push_back(r);

  int naxis = (iType==Cubic? 1: 2);
  for (int iaxis=0; iaxis<naxis; iaxis++)
  {
    cfgScalar Y = computeYoungModulus(strains[iaxis][iaxis], forceMagnitude);
    ioParameters.push_back(Y);
  }
  cfgScalar poissonRatio_xy = computePoissonRatio(strains[0][0], strains[0][1]); 
  ioParameters.push_back(poissonRatio_xy);

  cfgScalar G = C(2,2);
  ioParameters.push_back(G);

  delete stepper;
  delete physicalSystem;
}

void NumericalCoarsening::computeCoarsenedElasticityTensorAndParameters(const std::vector<int> &iMaterials, int N[2], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad2D> &iBaseMaterials, 
                                                                        StructureType iType, std::vector<cfgScalar> &iBaseMaterialDensities, std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters)
{
  std::vector<int> cellMaterials;
  repMaterialAssignment(N[0], N[1], iMaterials, iNbBlockRep, iNbBlockRep, iNbSubdiv, cellMaterials); 

  int n[2] = {N[0]*iNbBlockRep*iNbSubdiv, N[1]*iNbBlockRep*iNbSubdiv};
  ElementRegGrid2D * physicalSystem = createPhysicalSystem(n, iBaseMaterials, cellMaterials);
  Stepper2D * stepper = createStepper(physicalSystem);

  cfgScalar forceMagnitude = 1;
  std::vector<std::vector<cfgScalar> > harmonicDisp;
  computeHarmonicDisplacements(physicalSystem, stepper, forceMagnitude, iType, harmonicDisp);

  std::vector<Vector2S> h = cfgMaterialUtilities::toVector2S(harmonicDisp[0]);
  std::vector<Vector2S> x = cfgUtil::add(physicalSystem->X, h);

  MatrixXS C;
  meshUtil::computeCoarsenedElasticityTensor(*physicalSystem, harmonicDisp, C);

  std::vector<float> Cvalues;
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      if (j<=i)
      {
        ioTensorValues.push_back(C(i,j));
      }
    }
  }

  std::vector<std::vector<cfgScalar> > strains;
  meshUtil::computeStrains(*physicalSystem, harmonicDisp, strains);

  cfgScalar r = computeMaterialDensity(iMaterials, iBaseMaterialDensities);
  ioParameters.push_back(r);

  int naxis = (iType==Cubic? 1: 2);
  for (int iaxis=0; iaxis<naxis; iaxis++)
  {
    cfgScalar Y = computeYoungModulus(strains[iaxis][iaxis], forceMagnitude);
    ioParameters.push_back(Y);
  }
  cfgScalar poissonRatio_xy = computePoissonRatio(strains[0][0], strains[0][1]); 
  ioParameters.push_back(poissonRatio_xy);

  cfgScalar G = C(2,2);
  ioParameters.push_back(G);

  delete stepper;
  delete physicalSystem;
}

void NumericalCoarsening::computeCoarsenedElasticityTensorAndParameters3D(const std::vector<int> &iMaterials, int N[3], const std::vector<cfgScalar> &iBaseMaterialDensities, StructureType iType, std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters)
{
  std::vector<int> cellMaterials;
  repMaterialAssignment(N[0], N[1], N[2], iMaterials, m_nbBlockRep, m_nbBlockRep, m_nbBlockRep, m_nbSubdiv, cellMaterials); 

  int n[3] = {N[0]*m_nbBlockRep*m_nbSubdiv, N[1]*m_nbBlockRep*m_nbSubdiv, N[2]*m_nbBlockRep*m_nbSubdiv};
  ElementRegGrid * physicalSystem = createPhysicalSystem(n, m_baseMaterials3D, cellMaterials);

  int indVoidMaterial = std::min_element(iBaseMaterialDensities.begin(), iBaseMaterialDensities.end()) - iBaseMaterialDensities.begin();
  if (iBaseMaterialDensities[indVoidMaterial]==0 && m_baseMaterials3D[indVoidMaterial].e[0]->param[0]==0)
  {
    std::vector<int> isolatedPoints = getIsolatedPoints(cellMaterials, n[0], n[1], n[2], indVoidMaterial);
    setSubVector<int>(physicalSystem->fixed, 1, isolatedPoints);
  }

  cfgScalar forceMagnitude = 1;
  std::vector<std::vector<cfgScalar> > harmonicDisp;
  computeHarmonicDisplacements(physicalSystem, forceMagnitude, iType, harmonicDisp);

  MatrixXS C;
  meshUtil::computeCoarsenedElasticityTensor(*physicalSystem, harmonicDisp, C);
  //std::cout << "C = " << C << std::endl << std::endl;

  std::vector<float> Cvalues;
  for (int i=0; i<6; i++)
  {
    for (int j=0; j<6; j++)
    {
      if (j<=i)
      {
        ioTensorValues.push_back(C(i,j));
      }
    }
  }

  std::vector<std::vector<cfgScalar> > strains;
  meshUtil::computeStrains(*physicalSystem, harmonicDisp, strains);

  int nmat = (int)m_K0.size();
  std::vector<std::vector<int> > baseMaterials(nmat);
  for (int imat=0; imat<nmat; imat++)
  {
    baseMaterials[imat].push_back(imat);
  }

  //cfgScalar r = computeMaterialRatio(iMaterials, baseMaterials);
  cfgScalar r = computeMaterialDensity(iMaterials, iBaseMaterialDensities);
  ioParameters.push_back(r);

  int naxis = (iType==Cubic? 1: 3);
  for (int iaxis=0; iaxis<naxis; iaxis++)
  {
    cfgScalar Y = computeYoungModulus(strains[iaxis][iaxis], forceMagnitude);
    ioParameters.push_back(Y);
  }

  cfgScalar poissonRatio_xy = computePoissonRatio(strains[0][0], strains[0][1]); 
  ioParameters.push_back(poissonRatio_xy);

  if (iType!=Cubic)
  {
    cfgScalar poissonRatio_xz = computePoissonRatio(strains[0][0], strains[0][2]); 
    ioParameters.push_back(poissonRatio_xz);

    cfgScalar poissonRatio_yz = computePoissonRatio(strains[1][1], strains[1][2]); 
    ioParameters.push_back(poissonRatio_yz);
  }

  for (int i=0; i<naxis; i++)
  {
    cfgScalar G = C(3+i,3+i);
    ioParameters.push_back(G);
  }
  delete physicalSystem;
}

std::vector<int> NumericalCoarsening::getIsolatedPoints(const std::vector<int> &iMaterials, int nx, int ny, int indVoidMaterial)
{
  std::vector<int> isolatedPoints;

  for (int x=0; x<nx; x++)
  {
    for (int y=0; y<ny; y++)
    {
      int indMat11 = getGridToVectorIndex(x, y, nx, ny);
      int indMat12 = getGridToVectorIndex(x, (y+ny-1)%ny, nx, ny);
      int indMat21 = getGridToVectorIndex((x+nx-1)%nx, y, nx, ny);
      int indMat22 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, nx, ny);

      std::vector<int> materials;
      materials.push_back(iMaterials[indMat11]);
      materials.push_back(iMaterials[indMat12]);
      materials.push_back(iMaterials[indMat21]);
      materials.push_back(iMaterials[indMat22]);
      bool isolated = true;
      for (int imat=0; imat<(int)materials.size(); imat++)
      {
        int mat = materials[imat];
        isolated = isolated && (mat==indVoidMaterial);
      }
      if (isolated)
      {
        int indPoint = getGridToVectorIndex(x, y, nx+1, ny+1);
        isolatedPoints.push_back(indPoint);
      }
    }
  }
  return isolatedPoints;
}

std::vector<int> NumericalCoarsening::getIsolatedPoints(const std::vector<int> &iMaterials, int nx, int ny, int nz, int indVoidMaterial)
{
  std::vector<int> isolatedPoints;

  for (int x=0; x<nx; x++)
  {
    for (int y=0; y<ny; y++)
    {
      for (int z=0; z<nz; z++)
      {
        int indMat111 = getGridToVectorIndex(x, y, z, nx, ny, nz);
        int indMat121 = getGridToVectorIndex(x, (y+ny-1)%ny, z, nx, ny, nz);
        int indMat211 = getGridToVectorIndex((x+nx-1)%nx, y, z, nx, ny, nz);
        int indMat221 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, z, nx, ny, nz);
        int indMat112 = getGridToVectorIndex(x, y, (z+nz-1)%nz, nx, ny, nz);
        int indMat122 = getGridToVectorIndex(x, (y+ny-1)%ny, (z+nz-1)%nz, nx, ny, nz);
        int indMat212 = getGridToVectorIndex((x+nx-1)%nx, y, (z+nz-1)%nz, nx, ny, nz);
        int indMat222 = getGridToVectorIndex((x+nx-1)%nx, (y+ny-1)%ny, (z+nz-1)%nz, nx, ny, nz);

        std::vector<int> materials;
        materials.push_back(iMaterials[indMat111]);
        materials.push_back(iMaterials[indMat121]);
        materials.push_back(iMaterials[indMat211]);
        materials.push_back(iMaterials[indMat221]);
        materials.push_back(iMaterials[indMat112]);
        materials.push_back(iMaterials[indMat122]);
        materials.push_back(iMaterials[indMat212]);
        materials.push_back(iMaterials[indMat222]);
        bool isolated = true;
        for (int imat=0; imat<(int)materials.size(); imat++)
        {
          int mat = materials[imat];
          isolated = isolated && (mat==indVoidMaterial);
        }
        if (isolated)
        {
          int indPoint = getGridToVectorIndex(x, y, z, nx+1, ny+1, nz+1);
          isolatedPoints.push_back(indPoint);
        }
      }
    }
  }
  return isolatedPoints;
}


void NumericalCoarsening::computeCoarsenedElasticityTensorAndParameters2D(const std::vector<int> &iMaterials, int N[2], const std::vector<cfgScalar> &iBaseMaterialDensities, StructureType iType, std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters)
{
  std::vector<int> cellMaterials;
  repMaterialAssignment(N[0], N[1], iMaterials, m_nbBlockRep, m_nbBlockRep, m_nbSubdiv, cellMaterials); 

  int n[2] = {N[0]*m_nbBlockRep*m_nbSubdiv, N[1]*m_nbBlockRep*m_nbSubdiv};
  ElementRegGrid2D * physicalSystem = createPhysicalSystem(n, m_baseMaterials2D, cellMaterials);

  int indVoidMaterial = std::min_element(iBaseMaterialDensities.begin(), iBaseMaterialDensities.end()) - iBaseMaterialDensities.begin();
  if (iBaseMaterialDensities[indVoidMaterial]==0 && m_baseMaterials2D[indVoidMaterial].e[0]->param[0]==0)
  {
    std::vector<int> isolatedPoints = getIsolatedPoints(cellMaterials, n[0], n[1], indVoidMaterial);
    setSubVector<int>(physicalSystem->fixed, 1, isolatedPoints);
  }

  cfgScalar forceMagnitude = 1;
  std::vector<std::vector<cfgScalar> > harmonicDisp;
  computeHarmonicDisplacements(physicalSystem, forceMagnitude, iType, harmonicDisp);

  MatrixXS C;
  meshUtil::computeCoarsenedElasticityTensor(*physicalSystem, harmonicDisp, C);
  //std::cout << "C = " << C << std::endl << std::endl;

  std::vector<float> Cvalues;
  for (int i=0; i<3; i++)
  {
    for (int j=0; j<3; j++)
    {
      if (j<=i)
      {
        ioTensorValues.push_back(C(i,j));
      }
    }
  }

  std::vector<std::vector<cfgScalar> > strains;
  meshUtil::computeStrains(*physicalSystem, harmonicDisp, strains);

  int nmat = (int)m_K0.size();
  std::vector<std::vector<int> > baseMaterials(nmat);
  for (int imat=0; imat<nmat; imat++)
  {
    baseMaterials[imat].push_back(imat);
  }
  //cfgScalar r = computeMaterialRatio(iMaterials, baseMaterials);
  cfgScalar r = computeMaterialDensity(iMaterials, iBaseMaterialDensities);
  ioParameters.push_back(r);

  int naxis = (iType==Cubic? 1: 2);
  for (int iaxis=0; iaxis<naxis; iaxis++)
  {
    cfgScalar Y = computeYoungModulus(strains[iaxis][iaxis], forceMagnitude);
    ioParameters.push_back(Y);
  }

  cfgScalar poissonRatio_xy = computePoissonRatio(strains[0][0], strains[0][1]); 
  ioParameters.push_back(poissonRatio_xy);

  cfgScalar G = C(2,2);
  ioParameters.push_back(G);

  delete physicalSystem;
}

void NumericalCoarsening::computeCoarsenedElasticityTensorAndParameters(const std::vector<int> &iMaterials, int N[3], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad> &iBaseMaterials, 
                                                                        StructureType iType, std::vector<cfgScalar> &ioTensorValues, std::vector<cfgScalar> &ioParameters)
{
  std::vector<int> cellMaterials;
  repMaterialAssignment(N[0], N[1], N[1], iMaterials, iNbBlockRep, iNbBlockRep, iNbBlockRep, iNbSubdiv, cellMaterials); 

  int n[3] = {N[0]*iNbBlockRep*iNbSubdiv, N[1]*iNbBlockRep*iNbSubdiv, N[2]*iNbBlockRep*iNbSubdiv};
  ElementRegGrid * physicalSystem = createPhysicalSystem(n, iBaseMaterials, cellMaterials);
  Stepper * stepper = createStepper(physicalSystem);

  cfgScalar forceMagnitude = 1;
  std::vector<std::vector<cfgScalar> > harmonicDisp;
  computeHarmonicDisplacements(physicalSystem, stepper, forceMagnitude, iType, harmonicDisp);

  MatrixXS C;
  meshUtil::computeCoarsenedElasticityTensor(*physicalSystem, harmonicDisp, C);
  //std::cout << "C = " << C << std::endl << std::endl;

  std::vector<float> Cvalues;
  for (int i=0; i<6; i++)
  {
    for (int j=0; j<6; j++)
    {
      if (j<=i)
      {
        ioTensorValues.push_back(C(i,j));
      }
    }
  }

  std::vector<std::vector<cfgScalar> > strains;
  meshUtil::computeStrains(*physicalSystem, harmonicDisp, strains);

  std::vector<std::vector<int> > baseMaterials(2);
  baseMaterials[0].push_back(0);
  baseMaterials[1].push_back(1);

  cfgScalar r = computeMaterialRatio(iMaterials, baseMaterials);
  ioParameters.push_back(r);

  int naxis = (iType==Cubic? 1: 3);
  for (int iaxis=0; iaxis<naxis; iaxis++)
  {
    cfgScalar Y = computeYoungModulus(strains[iaxis][iaxis], forceMagnitude);
    ioParameters.push_back(Y);
  }

  cfgScalar poissonRatio_xy = computePoissonRatio(strains[0][0], strains[0][1]); 
  ioParameters.push_back(poissonRatio_xy);

  if (iType!=Cubic)
  {
    cfgScalar poissonRatio_xz = computePoissonRatio(strains[0][0], strains[0][2]); 
    ioParameters.push_back(poissonRatio_xz);

    cfgScalar poissonRatio_yz = computePoissonRatio(strains[1][1], strains[1][2]); 
    ioParameters.push_back(poissonRatio_yz);
  }

  for (int i=0; i<naxis; i++)
  {
    cfgScalar G = C(3+i,3+i);
    ioParameters.push_back(G);
  }
  delete stepper;
  delete physicalSystem;
}

MatrixXS NumericalCoarsening::computeCoarsenedElasticityTensor(const std::vector<int> &iMaterials, int N[2], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad2D> &iBaseMaterials)
{
  std::vector<int> cellMaterials;
  repMaterialAssignment(N[0], N[1], iMaterials, iNbBlockRep, iNbBlockRep, iNbSubdiv, cellMaterials); 

  int n[2] = {N[0]*iNbBlockRep*iNbSubdiv, N[1]*iNbBlockRep*iNbSubdiv};

  ElementRegGrid2D * physicalSystem = createPhysicalSystem(n, iBaseMaterials, cellMaterials);
  Stepper2D * stepper = createStepper(physicalSystem);

  cfgScalar forceMagnitude = 1;
  std::vector<std::vector<cfgScalar> > harmonicDisp;
  computeHarmonicDisplacements(physicalSystem, stepper, forceMagnitude, General, harmonicDisp);

  MatrixXS C;
  meshUtil::computeCoarsenedElasticityTensor(*physicalSystem, harmonicDisp, C);

  delete stepper;
  delete physicalSystem;

  return C;
}

MatrixXS NumericalCoarsening::computeCoarsenedElasticityTensor(const std::vector<int> &iMaterials, int N[3], int iNbBlockRep, int iNbSubdiv, std::vector<MaterialQuad> &iBaseMaterials)
{
  std::vector<int> cellMaterials;
  repMaterialAssignment(N[0], N[1], N[1], iMaterials, iNbBlockRep, iNbBlockRep, iNbBlockRep, iNbSubdiv, cellMaterials); 

  int n[3] = {N[0]*iNbBlockRep*iNbSubdiv, N[1]*iNbBlockRep*iNbSubdiv, N[2]*iNbBlockRep*iNbSubdiv};
  ElementRegGrid * physicalSystem = createPhysicalSystem(n, iBaseMaterials, cellMaterials);
  Stepper * stepper = createStepper(physicalSystem);
  
  cfgScalar forceMagnitude = 1;
  std::vector<std::vector<cfgScalar> > harmonicDisp;
  computeHarmonicDisplacements(physicalSystem, stepper, forceMagnitude, General, harmonicDisp);

  MatrixXS C;
  meshUtil::computeCoarsenedElasticityTensor(*physicalSystem, harmonicDisp, C);

  delete stepper;
  delete physicalSystem;

  return C;
}

ElementRegGrid2D * NumericalCoarsening::createPhysicalSystem(int n[2], std::vector<MaterialQuad2D> &iMaterials, const std::vector<int> &iMatAssignments)
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

ElementRegGrid * NumericalCoarsening::createPhysicalSystem(int n[3], std::vector<MaterialQuad> &iMaterials, const std::vector<int> &iMatAssignments)
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

Stepper2D * NumericalCoarsening::createStepper(ElementRegGrid2D * physicalSystem)
{
  Stepper2D * stepper = new StepperNewton2D();
  ((StepperNewton2D*)stepper)->removeTranslation(true);
  ((StepperNewton2D*)stepper)->removeRotation(true);
  ((StepperNewton2D*)stepper)->enforcePeriodicity(true);

  stepper->nSteps = 1;

  stepper->init(physicalSystem);

  return stepper;
}

Stepper* NumericalCoarsening::createStepper(ElementRegGrid * physicalSystem)
{
  Stepper* stepper = new StepperNewton();
  ((StepperNewton*)stepper)->removeTranslation(true);
  ((StepperNewton*)stepper)->removeRotation(true);
  ((StepperNewton*)stepper)->enforcePeriodicity(true);

  stepper->nSteps = 1;

  stepper->init(physicalSystem);

  return stepper;
}

void NumericalCoarsening::setExternalForcesForHarmonicDisp(cfgScalar iFx, cfgScalar iFy, ElementRegGrid2D *iPhysicalSystem)
{
  assert(iPhysicalSystem);

  int nforce =  iPhysicalSystem->fe.size();
  iPhysicalSystem->fe.clear();
  iPhysicalSystem->fe.resize(nforce, Vector2S(0,0));

  Vector2S force(iFx, iFy);
  if (iFx!=0 && iFy==0)
  {
    addExternalForces(iPhysicalSystem, 0, 2, force);
  }
  else if (iFx==0 && iFy!=0)
  {
    addExternalForces(iPhysicalSystem, 1, 2, force);
  }
  else if (iFx!=0 && iFy!=0)
  {
    Vector2S forceX(iFx, 0);
    addExternalForces(iPhysicalSystem, 1, 2, forceX);

    Vector2S forceY(0, iFy);
    addExternalForces(iPhysicalSystem, 0, 2, forceY);
  }
}

void NumericalCoarsening::setExternalForcesForHarmonicDisp(cfgScalar iFx, cfgScalar iFy, cfgScalar iFz, ElementRegGrid *iPhysicalSystem)
{
  assert(iPhysicalSystem);

  int nforce =  iPhysicalSystem->fe.size();
  iPhysicalSystem->fe.clear();
  iPhysicalSystem->fe.resize(nforce, Vector3S(0,0,0));

  Vector3S force(iFx, iFy, iFz);
  if (iFx!=0 && iFy==0 && iFz==0)
  {
    addExternalForces(iPhysicalSystem, 0, 2, force);
  }
  else if (iFx==0 && iFy!=0 && iFz==0)
  {
    addExternalForces(iPhysicalSystem, 1, 2, force);
  }
  else if (iFx==0 && iFy==0 && iFz!=0)
  {
    addExternalForces(iPhysicalSystem, 2, 2, force);
  }
  else if (iFx!=0 && iFy!=0 && iFz==0)
  {
    Vector3S forceX(iFx, 0, 0);
    addExternalForces(iPhysicalSystem, 1, 2, forceX);

    Vector3S forceY(0, iFy, 0);
    addExternalForces(iPhysicalSystem, 0, 2, forceY);
  }
  else if (iFx!=0 && iFy==0 && iFz!=0)
  {
    Vector3S forceX(iFx, 0, 0);
    addExternalForces(iPhysicalSystem, 2, 2, forceX);

    Vector3S forceZ(0, 0, iFz);
    addExternalForces(iPhysicalSystem, 0, 2, forceZ);
  }
  else if (iFx==0 && iFy!=0 && iFz!=0)
  {
    Vector3S forceY(0, iFy, 0);
    addExternalForces(iPhysicalSystem, 2, 2, forceY);

    Vector3S forceZ(0, 0, iFz);
    addExternalForces(iPhysicalSystem, 1, 2, forceZ);
  }
}

void NumericalCoarsening::getExternalForcesForHarmonicDisp(cfgScalar iFx, cfgScalar iFy, ElementRegGrid2D *iPhysicalSystem, std::vector<float> &oForces)
{
  assert(iPhysicalSystem);

  Vector2S force(iFx, iFy);
  if (iFx!=0 && iFy==0)
  {
    getExternalForces(iPhysicalSystem, 0, 2, force, oForces);
  }
  else if (iFx==0 && iFy!=0)
  {
    getExternalForces(iPhysicalSystem, 1, 2, force, oForces);
  }
  else if (iFx!=0 && iFy!=0)
  {
    Vector2S forceX(iFx, 0);
    getExternalForces(iPhysicalSystem, 1, 2, forceX, oForces);

    Vector2S forceY(0, iFy);
    getExternalForces(iPhysicalSystem, 0, 2, forceY, oForces);
  }
}

void NumericalCoarsening::getExternalForcesForHarmonicDisp(cfgScalar iFx, cfgScalar iFy, cfgScalar iFz, ElementRegGrid *iPhysicalSystem, std::vector<float> &oForces)
{
  assert(iPhysicalSystem);

  Vector3S force(iFx, iFy, iFz);
  if (iFx!=0 && iFy==0 && iFz==0)
  {
    getExternalForces(iPhysicalSystem, 0, 2, force, oForces);
  }
  else if (iFx==0 && iFy!=0 && iFz==0)
  {
    getExternalForces(iPhysicalSystem, 1, 2, force, oForces);
  }
  else if (iFx==0 && iFy==0 && iFz!=0)
  {
    getExternalForces(iPhysicalSystem, 2, 2, force, oForces);
  }
  else if (iFx!=0 && iFy!=0 && iFz==0)
  {
    Vector3S forceX(iFx, 0, 0);
    getExternalForces(iPhysicalSystem, 1, 2, forceX, oForces);

    Vector3S forceY(0, iFy, 0);
    getExternalForces(iPhysicalSystem, 0, 2, forceY, oForces);
  }
  else if (iFx!=0 && iFy==0 && iFz!=0)
  {
    Vector3S forceX(iFx, 0, 0);
    getExternalForces(iPhysicalSystem, 2, 2, forceX, oForces);

    Vector3S forceZ(0, 0, iFz);
    getExternalForces(iPhysicalSystem, 0, 2, forceZ, oForces);
  }
  else if (iFx==0 && iFy!=0 && iFz!=0)
  {
    Vector3S forceY(0, iFy, 0);
    getExternalForces(iPhysicalSystem, 2, 2, forceY, oForces);

    Vector3S forceZ(0, 0, iFz);
    getExternalForces(iPhysicalSystem, 1, 2, forceZ, oForces);
  }
}

void NumericalCoarsening::addExternalForces(ElementRegGrid2D * iElementGrid, int iAxis, int iSide, Vector2S &iForceMagnitude)
{
  assert(iElementGrid);

  ElementRegGrid2D * em = iElementGrid;

  int n[2];
  n[0] = em->nx;
  n[1] = em->ny;

  cfgScalar thickness = 1.;
  cfgScalar coeff = ((cfgScalar)1/(n[(iAxis+1)%2]*thickness))/2;

  Vector2S ff = coeff*iForceMagnitude;

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
        em->fe[VertexIndex] += (cfgScalar)sign*ff;
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

void NumericalCoarsening::getExternalForces(ElementRegGrid2D * iElementGrid, int iAxis, int iSide, Vector2S &iForceMagnitude, std::vector<cfgScalar> &oForces)
{
  assert(iElementGrid);

  ElementRegGrid2D * em = iElementGrid;
  oForces.clear();
  oForces.resize(2*em->fe.size());

  int n[2];
  n[0] = em->nx;
  n[1] = em->ny;

  cfgScalar thickness = 1.;
  cfgScalar coeff = ((cfgScalar)1/(n[(iAxis+1)%2]*thickness))/2;

  Vector2S ff = coeff*iForceMagnitude;

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
        oForces[2*VertexIndex] += (cfgScalar)sign*ff[0];
        oForces[2*VertexIndex+1] += (cfgScalar)sign*ff[1];
      }
    }
  }
}

void NumericalCoarsening::getExternalForces(ElementRegGrid * iElementGrid, int iAxis, int iSide, Vector3S &iForceMagnitude, std::vector<cfgScalar> &oForces)
{
  assert(iElementGrid);
  assert(iSide==2);

  ElementRegGrid * em = iElementGrid;
  oForces.resize(3*em->fe.size());

  int n[3];
  n[0] = em->nx;
  n[1] = em->ny;
  n[2] = em->nz;

  std::vector<int> sides;
  std::vector<int> signs;

  cfgScalar coeff = ((cfgScalar)1/(n[(iAxis+1)%3]*n[(iAxis+2)%3]))/4;
  Vector3S ff = coeff*iForceMagnitude;

  sides.push_back(0);
  sides.push_back(1);

  signs.push_back(-1);
  signs.push_back(1);

  for (int iside=0; iside<2; iside++)
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
        oForces[3*VertexIndex] += (cfgScalar)sign*ff[0];
        oForces[3*VertexIndex+1] += (cfgScalar)sign*ff[1];
        oForces[3*VertexIndex+2] += (cfgScalar)sign*ff[2];
      }
    }
  }
}

void NumericalCoarsening::addExternalForces(ElementRegGrid * iElementGrid, int iAxis, int iSide, Vector3S &iForceMagnitude)
{
  assert(iElementGrid);
  assert(iSide==2);

  ElementRegGrid * em = iElementGrid;

  int n[3];
  n[0] = em->nx;
  n[1] = em->ny;
  n[2] = em->nz;

  std::vector<int> sides;
  std::vector<int> signs;

  cfgScalar coeff = ((cfgScalar)1/(n[(iAxis+1)%3]*n[(iAxis+2)%3]))/4;
  Vector3S ff = coeff*iForceMagnitude;

  sides.push_back(0);
  sides.push_back(1);

  signs.push_back(-1);
  signs.push_back(1);

  for (int iside=0; iside<2; iside++)
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
        em->fe[VertexIndex] += (cfgScalar)sign*ff;
      }
    }
  }
}

void NumericalCoarsening::computeHarmonicDisplacements(ElementRegGrid2D * iPhysicalSystem, Stepper2D * iStepper, cfgScalar iForceMagnitude, StructureType iType, std::vector<std::vector<cfgScalar> > &oHarmonicDisplacements)
{
  oHarmonicDisplacements.clear();

  cfgScalar forces[3][2] = { {iForceMagnitude,0}, {0,iForceMagnitude}, {iForceMagnitude/2, iForceMagnitude/2 }};

  for (int idisp=0; idisp<3; idisp++)
  {
    iPhysicalSystem->x = iPhysicalSystem->X;
    setExternalForcesForHarmonicDisp(forces[idisp][0], forces[idisp][1], iPhysicalSystem);
    iStepper->step();

    std::vector<cfgScalar> X = toVectorScalar(iPhysicalSystem->X);
    std::vector<cfgScalar> x = toVectorScalar(iPhysicalSystem->x);

    oHarmonicDisplacements.push_back(sub(x, X));
  }
  /*std::vector<std::vector<double> > setForces(3);
  for (int idisp=0; idisp<3; idisp++)
  {
    std::vector<cfgScalar> currentForces;
    getExternalForcesForHarmonicDisp(forces[idisp][0], forces[idisp][1], iPhysicalSystem, currentForces);
    setForces[idisp] = convertVec<cfgScalar, double>(currentForces);
  }*/ 
}

void NumericalCoarsening::computeHarmonicDisplacements(ElementRegGrid * iPhysicalSystem, Stepper * iStepper, cfgScalar iForceMagnitude, StructureType iType, std::vector<std::vector<cfgScalar> > &oHarmonicDisplacements)
{
  oHarmonicDisplacements.clear();

  cfgScalar forces[6][3] = { {iForceMagnitude, 0, 0}, {0, iForceMagnitude, 0}, {0, 0, iForceMagnitude}, {iForceMagnitude/2, iForceMagnitude/2, 0}, {0, iForceMagnitude/2, iForceMagnitude/2}, {iForceMagnitude/2, 0, iForceMagnitude/2}};

  /*
  for (int idisp=0; idisp<6; idisp++)
  {
    iPhysicalSystem->x = iPhysicalSystem->X;
    setExternalForcesForHarmonicDisp(forces[idisp][0], forces[idisp][1], forces[idisp][2], iPhysicalSystem);

    iStepper->step();

    std::vector<cfgScalar> X = toVectorFloat(iPhysicalSystem->X);
    std::vector<cfgScalar> x = toVectorFloat(iPhysicalSystem->x);

    oHarmonicDisplacements.push_back(sub(x, X));
  }*/ 

  std::vector<std::vector<double> > externalForces(6);
  for (int idisp=0; idisp<6; idisp++)
  {
    std::vector<cfgScalar> currentForces;
    getExternalForcesForHarmonicDisp(forces[idisp][0], forces[idisp][1], forces[idisp][2], iPhysicalSystem, currentForces);
    externalForces[idisp] = convertVec<cfgScalar, double>(currentForces);
  } 
  int nforce = externalForces[0].size();

  Timer t;
  std::vector<int> m_I, m_J;
  bool triangular = true;
  m_I.clear();
  m_J.clear();

  //t.start();
  iPhysicalSystem->stiffnessPattern(m_I, m_J, triangular, true, true, true);
  //t.end();
  //std::cout << " stiffnessPattern: " <<  t.getSeconds() << std::endl;
  //pardiso uses 1-based
  for (unsigned int ii = 0; ii < m_I.size(); ii++){
    m_I[ii] ++;
  }
  for (unsigned int ii = 0; ii < m_J.size(); ii++){
    m_J[ii] ++;
  }
  PardisoState * pardisoState = new PardisoState();
  //t.start();
  pardisoInit(pardisoState);
  //t.end();
  //std::cout << " pardisoInit: " <<  t.getSeconds() << std::endl;
  //t.start();
  pardisoSymbolicFactorize(m_I.data(), m_J.data(), (int)m_I.size()-1, pardisoState);
  //t.end();
  //std::cout << " pardisoSymbolicFactorize: " <<  t.getSeconds() << std::endl;

  //t.start();
  std::vector<float> Kvalues;
  iPhysicalSystem->getStiffnessSparse(Kvalues, triangular, true, true, true, true);
  std::vector<double> val = convertVec<cfgScalar, double>(Kvalues);
  //t.end();
  //std::cout << " getStiffnessSparse: " <<  t.getSeconds() << std::endl;

  //t.start();
  int nrows = (int)m_I.size() - 1;
  pardisoNumericalFactorize(m_I.data(), m_J.data(), val.data(), nrows, pardisoState);
  //t.end();
  //std::cout << " pardisoNumericalFactorize: " <<  t.getSeconds() << std::endl;

  for (unsigned int ii = 0; ii < externalForces.size(); ii++)
  {
    externalForces[ii].resize(nrows);
    std::vector<double> u(nrows, 0.);
    //t.start();
    pardisoBackSubstitute(m_I.data(), m_J.data(), val.data(), nrows, u.data(), externalForces[ii].data(), pardisoState);
    //t.end();
    //std::cout << " pardisoBackSubstitute: " <<  t.getSeconds() << std::endl;
    u.resize(nforce);
    oHarmonicDisplacements.push_back(convertVec<double, cfgScalar>(u));
  }
  pardisoFree(m_I.data(), m_J.data(), nrows, pardisoState);
  delete pardisoState; 
}

void NumericalCoarsening::computeHarmonicDisplacements(ElementRegGrid * iPhysicalSystem, cfgScalar iForceMagnitude, StructureType iType, std::vector<std::vector<cfgScalar> > &oHarmonicDisplacements)
{
  oHarmonicDisplacements.clear();

  cfgScalar forces[6][3] = { {iForceMagnitude, 0, 0}, {0, iForceMagnitude, 0}, {0, 0, iForceMagnitude}, {iForceMagnitude/2, iForceMagnitude/2, 0}, {0, iForceMagnitude/2, iForceMagnitude/2}, {iForceMagnitude/2, 0, iForceMagnitude/2}};

  std::vector<std::vector<double> > externalForces(6);
  for (int idisp=0; idisp<6; idisp++)
  {
    std::vector<cfgScalar> currentForces;
    getExternalForcesForHarmonicDisp(forces[idisp][0], forces[idisp][1], forces[idisp][2], iPhysicalSystem, currentForces);
    externalForces[idisp] = convertVec<cfgScalar, double>(currentForces);
  } 
  int nforce = externalForces[0].size();

  bool triangular = true;
  //Timer t;
  //t.start();
  std::vector<double> val;
  getStiffnessSparse(iPhysicalSystem, val);
  std::vector<float> Kvalues;
  //iPhysicalSystem->getStiffnessSparse(Kvalues, triangular, true, true, true, true);
  //std::vector<double> val = convertVec<cfgScalar, double>(Kvalues);
  //t.end();
  //std::cout << " getStiffnessSparse: " <<  t.getSeconds() << std::endl;

  //t.start();
  int nrows = (int)m_I.size() - 1;
  pardisoNumericalFactorize(m_I.data(), m_J.data(), val.data(), nrows, m_pardisoState);
  //t.end();
  //std::cout << " pardisoNumericalFactorize: " <<  t.getSeconds() << std::endl;

  for (unsigned int ii = 0; ii < externalForces.size(); ii++)
  {
    externalForces[ii].resize(nrows);
    std::vector<double> u(nrows, 0.);
    //t.start();
    pardisoBackSubstitute(m_I.data(), m_J.data(), val.data(), nrows, u.data(), externalForces[ii].data(), m_pardisoState);
    //t.end();
    //std::cout << " pardisoBackSubstitute: " <<  t.getSeconds() << std::endl;
    u.resize(nforce);
    oHarmonicDisplacements.push_back(convertVec<double, cfgScalar>(u));
  }
}

void NumericalCoarsening::computeHarmonicDisplacements(ElementRegGrid2D * iPhysicalSystem, cfgScalar iForceMagnitude, StructureType iType, std::vector<std::vector<cfgScalar> > &oHarmonicDisplacements)
{
  oHarmonicDisplacements.clear();

  cfgScalar forces[3][2] = { {iForceMagnitude,0}, {0,iForceMagnitude}, {iForceMagnitude/2, iForceMagnitude/2 }};

  std::vector<std::vector<double> > externalForces(3);
  for (int idisp=0; idisp<3; idisp++)
  {
    std::vector<cfgScalar> currentForces;
    getExternalForcesForHarmonicDisp(forces[idisp][0], forces[idisp][1], iPhysicalSystem, currentForces);
    externalForces[idisp] = convertVec<cfgScalar, double>(currentForces);
  } 
  int nforce = externalForces[0].size();

  bool triangular = true;
  //Timer t;
  //t.start();
  std::vector<double> val;
  getStiffnessSparse(iPhysicalSystem, val);
  std::vector<float> Kvalues;
  //iPhysicalSystem->getStiffnessSparse(Kvalues, triangular, true, true, true, true);
  //std::vector<double> val = convertVec<cfgScalar, double>(Kvalues);
  //t.end();
  //std::cout << " getStiffnessSparse: " <<  t.getSeconds() << std::endl;

  //t.start();
  int nrows = (int)m_I.size() - 1;
  pardisoNumericalFactorize(m_I.data(), m_J.data(), val.data(), nrows, m_pardisoState);
  //t.end();
  //std::cout << " pardisoNumericalFactorize: " <<  t.getSeconds() << std::endl;

  for (unsigned int ii = 0; ii < externalForces.size(); ii++)
  {
    externalForces[ii].resize(nrows);
    std::vector<double> u(nrows, 0.);
    //t.start();
    pardisoBackSubstitute(m_I.data(), m_J.data(), val.data(), nrows, u.data(), externalForces[ii].data(), m_pardisoState);
    //t.end();
    //std::cout << " pardisoBackSubstitute: " <<  t.getSeconds() << std::endl;
    u.resize(nforce);
    oHarmonicDisplacements.push_back(convertVec<double, cfgScalar>(u));
  }
}

void NumericalCoarsening::getStiffnessSparse(ElementRegGrid * iPhysicalSystem, std::vector<double> &oValues)
{
  bool fixRigid = true;
  bool periodic = true;
  bool constrained = true;
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
            if (constrained)
            {
              if (iPhysicalSystem->fixed[vk] || iPhysicalSystem->fixed[vj])
              {
                if (vj != vk || dim1 != dim2)
                {
                  val = 0;
                }
                else
                {
                  val = 1;
                }
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

void NumericalCoarsening::getStiffnessSparse(ElementRegGrid2D * iPhysicalSystem, std::vector<double> &oValues)
{
  bool fixRigid = true;
  bool periodic = true;
  bool constrained = true;
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
            if (constrained)
            {
              if (iPhysicalSystem->fixed[vk] || iPhysicalSystem->fixed[vj]) 
              {
                if (vj != vk || dim1 != dim2)
                {
                  val = 0;
                }
                else
                {
                  val = 1;
                }
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
