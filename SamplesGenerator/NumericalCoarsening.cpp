
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

#include "MeshUtilities.h"
using namespace meshUtil;

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

#include "cfgUtilities.h"
using namespace cfgUtil;

NumericalCoarsening::NumericalCoarsening()
{
}

NumericalCoarsening::~NumericalCoarsening()
{
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
}

void NumericalCoarsening::computeHarmonicDisplacements(ElementRegGrid * iPhysicalSystem, Stepper * iStepper, cfgScalar iForceMagnitude, StructureType iType, std::vector<std::vector<cfgScalar> > &oHarmonicDisplacements)
{
  oHarmonicDisplacements.clear();

  cfgScalar forces[6][3] = { {iForceMagnitude, 0, 0}, {0, iForceMagnitude, 0}, {0, 0, iForceMagnitude}, {iForceMagnitude/2, iForceMagnitude/2, 0}, {0, iForceMagnitude/2, iForceMagnitude/2}, {iForceMagnitude/2, 0, iForceMagnitude/2}};

  for (int idisp=0; idisp<6; idisp++)
  {
    iPhysicalSystem->x = iPhysicalSystem->X;
    setExternalForcesForHarmonicDisp(forces[idisp][0], forces[idisp][1], forces[idisp][2], iPhysicalSystem);

    iStepper->step();

    std::vector<cfgScalar> X = toVectorFloat(iPhysicalSystem->X);
    std::vector<cfgScalar> x = toVectorFloat(iPhysicalSystem->x);

    oHarmonicDisplacements.push_back(sub(x, X));
  }
}

