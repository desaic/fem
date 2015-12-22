#include "MaterialOptimizer.h"

#include "EigenUtil.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "ElementRegGrid2D.h"
#include "Element2D.h"
#include "FEM2DFun.hpp"
#include "PiecewiseConstant2D.hpp"
#include "PiecewiseConstantSym2D.hpp"

#include "MaterialOptimizer.h"

#include "MaterialQuad2D.h"
#include "StrainLin2D.h"

#include "MeshUtilities.h"
using namespace meshUtil;

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

MaterialOptimizer::MaterialOptimizer()
{
  m_structureType = General;
}

MaterialOptimizer::~MaterialOptimizer()
{

}

void MaterialOptimizer::setBaseMaterials(const std::vector<MaterialQuad2D> &iBaseMaterials)
{
  m_baseMaterials = iBaseMaterials;
}

void MaterialOptimizer::setStructureType(StructureType iType)
{
  m_structureType = iType;
}

bool MaterialOptimizer::run(int N[2], std::vector<int> &ioNewMaterials, const std::vector<float> &iTargetParams)
{ 
  bool resOk = true;

  int nx = N[0], ny = N[1];

  ElementRegGrid2D * em = new ElementRegGrid2D(nx,ny);
  std::vector<StrainLin2D> ene(1);

  // Materials
  /*int imat=0, nmat=(int)m_baseMaterials.size();
  for (imat=0; imat<nmat; imat++)
  {
    em->addMaterial(&m_baseMaterials[imat]);
  }*/ 
  em->addMaterial(&m_baseMaterials[1]);
  em->initArrays();
  em->check();

  //assert(em->me.size()==ioNewMaterials.size());
  //em->me = ioNewMaterials;

  FEM2DFun * fem = new FEM2DFun();
  RealField * field = NULL;
  if (m_structureType == Cubic)
  {
    //PiecewiseConstant2D * field = new PiecewiseConstant2D();
    PiecewiseConstantSym2D * symfield = new PiecewiseConstantSym2D();
    symfield->allocate(nx/2, ny/2);
    field = symfield;
  }
  fem->field = field;
  fem->lowerBounds = 1e-3 * Eigen::VectorXd::Ones(field->param.size());
  fem->upperBounds = Eigen::VectorXd::Ones(field->param.size());

  fem->em = em;
  fem->m_nx = nx;
  fem->m_ny = ny;

  assert(iTargetParams.size()==4);
  float Y = iTargetParams[1];
  float nu = iTargetParams[2];

  float fx = 1;
  float dx0, dy0;
  computeTargetDisplacements(Y, nu, fx, dx0, dy0);

  fem->dx0 = (double)dx0;
  fem->dy0 = (double)dy0;
  fem->dxw = 1;
  fem->dyw = 1;

  std::vector< std::vector<double> > externalForces(1);
  getExternalForces(em, 0, 2, Vector2S(fx, 0), externalForces[0]);
  fem->externalForce = externalForces;
  fem->forceMagnitude = (double)fx;

  std::vector<int> matAssignment;
  if (m_structureType == Cubic)
  {
    getQuarter(nx, ny, ioNewMaterials, matAssignment);
  }
  else
  {
    matAssignment = ioNewMaterials;
  }

  double shrink = 0.1;
  Eigen::VectorXd x0(field->param.size());
  assert(x0.rows()==(int)matAssignment.size());
  for (int i=0; i<x0.rows(); i++)
  {
    //x0[i] = (float)matAssignment[i];
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
  fem->init(x0);

  double h = 1e-2;
  //check_df(fem, x0, h);
  fem->setParam(x0);
  Eigen::VectorXd grad = fem->df();
  h = 1;
  gradientDescent(fem, x0, 1000);

  for (int i=0; i<x0.rows(); i++)
  {
    if (x0[i] < 0.5)
    {
      matAssignment[i] = 0;
    }
    else
    {
      matAssignment[i] = 1;
    }
  }
  if (m_structureType == Cubic)
  {
    mirrorStructure(nx/2, ny/2, matAssignment, ioNewMaterials);
  }
  else
  {
    ioNewMaterials = matAssignment;
  }
  return resOk;
}

double MaterialOptimizer::infNorm(const Eigen::VectorXd & a)
{
  double maxval = 0;
  for (int ii = 0; ii < a.rows(); ii++){
    maxval = std::max(maxval, std::abs(a[ii]));
  }
  return maxval;
}

void MaterialOptimizer::gradientDescent(RealFun * fun, Eigen::VectorXd & x0, int nSteps)
{
  //maximum movement in any parameter.
  double maxStep = 1.;
  Eigen::VectorXd x = x0;
  for (int ii = 0; ii < nSteps; ii++){
    fun->setParam(x);
    Eigen::VectorXd grad = fun->df();
    double h = 1;
    double norm = infNorm(grad);
    h = maxStep / norm;
    int ret = lineSearch(fun, x, grad, h);
    //std::cout << ii << " " << fun->f() << " " << h << "\n";
    if (ret < 0){
      break;
    }
    x0 = x;
  }
}

int MaterialOptimizer::lineSearch(RealFun * fun, Eigen::VectorXd & x0, const Eigen::VectorXd & dir, double & h )
{
  double f0 = fun->f();
  Eigen::VectorXd grad0 = fun->df();
  double norm0 = infNorm(grad0);
  while (h>1.e-8){
    //minus sign here if we want to minimize function value.
    Eigen::VectorXd x = x0 - h*dir;
    clampVector(x, fun->lowerBounds, fun->upperBounds);
    fun->setParam(x);
    double f1 = fun->f();
    Eigen::VectorXd grad1 = fun->df();
    double norm1 = infNorm(grad1);
    //if gradient norm or function value decrease, return.
    if (norm1 < norm0 || f1 < f0){
      x0 = x;
      break;
    }
    else{
      h /= 2;
    }
  }
  return 0;
}

void MaterialOptimizer::check_df(RealFun * fun, const Eigen::VectorXd & x0, double h)
{
  fun->setParam(x0);
  Eigen::VectorXd x = x0;
  Eigen::VectorXd ana_df = fun->df();
  Eigen::VectorXd num_df = Eigen::VectorXd::Zero(x0.size());

  double max_diff = 0;
  for (int ii = 0; ii < x0.rows(); ii++){
    x[ii] += h;
    fun->setParam(x);
    double f_plus = fun->f();

    x[ii] -= 2 * h;
    fun->setParam(x);
    double f_minus = fun->f();

    x[ii] += h;

    num_df[ii] = (f_plus - f_minus) / (2 * h);
    std::cout << ii << " " << ana_df[ii] << " " << num_df[ii] << "\n";
    double diff = ana_df[ii] - num_df[ii];
    max_diff = std::max(max_diff, std::abs(diff));
  }
  std::cout << "max diff " << max_diff << "\n";
}

void MaterialOptimizer::computeTargetDisplacements(float iYoungModulus, float iPoissonsRatio, float iFx, float &odx, float &ody)
{
  //Y = f/dx
  odx = iFx/iYoungModulus;

  //nu = -dy/dx
  ody = -odx*iPoissonsRatio;
}

void MaterialOptimizer::getExternalForces(ElementRegGrid2D * iElementGrid, int iAxis, int iSide, Vector2S &iForceMagnitude, std::vector<double> &oForces)
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
        int VertexIndex = em->e[indElement]->at(FvIndex);
        oForces[2*VertexIndex] += (double)sign*ff[0];
        oForces[2*VertexIndex+1] += (double)sign*ff[1];
      }
    }
  }
}


