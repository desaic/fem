#include "MaterialOptimizer.h"

#include "EigenUtil.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "ElementRegGrid2D.h"
#include "Element2D.h"
#include "FEM2DFun.hpp"
#include "FEM3DFun.hpp"
#include "PiecewiseConstant2D.hpp"
#include "PiecewiseConstantSym2D.hpp"
#include "PiecewiseConstantCubic2D.hpp"
#include "PiecewiseConstant3D.hpp"
#include "PiecewiseConstantSym3D.hpp"
#include "PiecewiseConstantCubic3D.hpp"

#include "MaterialOptimizer.h"
#include <set>

#include "MaterialQuad2D.h"
#include "MaterialQuad.hpp"
#include "StrainLin2D.h"
#include "StrainLin.hpp"

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

void MaterialOptimizer::setStructureType(StructureType iType)
{
  m_structureType = iType;
}

bool MaterialOptimizer::run2D(int N[2], const std::vector<MaterialQuad2D> &iBaseMaterials, const std::vector<int> &iMaterialAssignments, const std::vector<float> &iTargetParams, std::vector<std::vector<int> > &oNewMaterialAssignments)
{ 
  oNewMaterialAssignments.clear();

  bool resOk = true;

  int nx = N[0], ny = N[1];

  ElementRegGrid2D * em = new ElementRegGrid2D(nx,ny);
  std::vector<StrainLin2D> ene(1);

  std::vector<MaterialQuad2D> baseMaterials = iBaseMaterials;
  em->addMaterial(&baseMaterials[1]);
  em->initArrays();
  em->check();

  FEM2DFun * fem = new FEM2DFun();
  RealField * field = NULL;
  if (m_structureType == Cubic)
  {
    PiecewiseConstantCubic2D * cubicfield = new PiecewiseConstantCubic2D();
    cubicfield->allocate(nx/2, ny/2);
    field = cubicfield;
  }
  else if (m_structureType == Orthotropic)
  {
    PiecewiseConstantSym2D * symfield = new PiecewiseConstantSym2D();
    symfield->allocate(nx/2, ny/2);
    field = symfield;
  }
  fem->field = field;
  fem->lowerBounds = 1e-3 * Eigen::VectorXd::Ones(field->param.size());
  fem->upperBounds = Eigen::VectorXd::Ones(field->param.size());

  fem->em = em;
  fem->gridSize[0] = nx;
  fem->gridSize[1] = ny;

  std::vector<int> matAssignment;
  if (m_structureType == Cubic || m_structureType == Orthotropic)
  {
    getQuarter(nx, ny, iMaterialAssignments, matAssignment);
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
  fem->setParam(x0);

  Eigen::MatrixXd target_G ;
  if (m_structureType == Cubic)
  {
    assert(iTargetParams.size()==4);
    float E = iTargetParams[1];
    float nu = iTargetParams[2];
    float mu = (float)fem->G(2,2);
    target_G = computeTargetStrains(2, E, nu, mu);
  }
  else if (m_structureType == Orthotropic)
  {
    assert(iTargetParams.size()==5);
    float E_x = iTargetParams[1];
    float E_y = iTargetParams[2];
    float nu_xy = iTargetParams[3];
    float mu_xy = (float)fem->G(2,2);
    target_G = computeTargetStrains(E_x, E_y, nu_xy, mu_xy);
  }
  else
  {
    // not implemented
    assert(0);
  }
  std::cout << "target G = " << std::endl;
  std::cout << target_G << std::endl << std::endl;

  float fx = 1;
  fem->forceMagnitude = (double)fx;

  fem->m0 = 0.5 * sum(fem->distribution) / fem->distribution.size();
  fem->mw = 0;//0.1 * fem->G(0, 0) / fem->density;
  double weightMu = 0.;
  fem->G0 = target_G;
  fem->wG(0) = 1;
  fem->wG(1) = 1;
  fem->wG(2) = weightMu;

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
        /*std::vector<int> newMatAssignment;
        mirrorStructure(nx/2, ny/2, mat, newMatAssignment);
        oNewMaterialAssignments.push_back(newMatAssignment);*/ 
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

bool MaterialOptimizer::run3D(int N[3], const std::vector<MaterialQuad> &iBaseMaterials, const std::vector<int> &iMaterialAssignments, const std::vector<float> &iTargetParams, std::vector<std::vector<int> > &oNewMaterialAssignments)
{ 
  oNewMaterialAssignments.clear();

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
  fem->init(x0);
  fem->setParam(x0);
  double val = fem->f();

  assert(iTargetParams.size()==4);
  float E = iTargetParams[1];
  float nu = iTargetParams[2];
  float mu = (float)fem->G(3,3);

  Eigen::MatrixXd target_G = computeTargetStrains(3, E, nu, mu);
  std::cout << "target G = " << std::endl;
  std::cout << target_G << std::endl << std::endl;

  float fx = 1;
  fem->forceMagnitude = (double)fx;

  fem->m0 = 0.5 * sum(fem->distribution) / fem->distribution.size();
  fem->mw = 0; //0.1 * fem->G(0, 0) / fem->density;
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
        /*std::vector<int> newMatAssignment;
        mirrorStructure(nx/2, ny/2, mat, newMatAssignment);
        oNewMaterialAssignments.push_back(newMatAssignment);*/ 
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

double MaterialOptimizer::infNorm(const Eigen::VectorXd & a)
{
  double maxval = 0;
  for (int ii = 0; ii < a.rows(); ii++){
    maxval = std::max(maxval, std::abs(a[ii]));
  }
  return maxval;
}

void MaterialOptimizer::gradientDescent(FEM2DFun * fun, Eigen::VectorXd & x0, int nSteps, std::vector<Eigen::VectorXd> &oParameters)
{
  //oParameters.push_back(x0);
  oParameters.push_back(fun->distribution);
  //maximum movement in any parameter.
  double maxStep = 0.5;
  Eigen::VectorXd x = x0;
  double h = 1;
  for (int ii = 0; ii < nSteps; ii++)
  {
    if (ii%50==0)
    {
      std::cout << "step = " << ii << std::endl;
    }
    fun->setParam(x);
    double val = fun->f();
    Eigen::VectorXd grad = fun->df();
    double norm = infNorm(grad);
    if (ii==0)
    {
      h = maxStep / norm;
    }
    double init_h = h;
    int ret = lineSearch(fun, x, grad, h);
    std::cout << "step = " << ii << "                              " << " f = " << val << " norm = " << norm << " " << "h = " << init_h << " " << h << std::endl;
    if (ret < 0){
      break;
    }
    x0 = x;
    //oParameters.push_back(x0);
    oParameters.push_back(fun->distribution);
  }
}

void MaterialOptimizer::gradientDescent(FEM3DFun * fun, Eigen::VectorXd & x0, int nSteps, std::vector<Eigen::VectorXd> &oParameters)
{
  //oParameters.push_back(x0);
  oParameters.push_back(fun->distribution);
  //maximum movement in any parameter.
  double maxStep = 0.5;
  Eigen::VectorXd x = x0;
  double h = 1;

  for (int ii = 0; ii < nSteps; ii++)
  {
    if (ii%50==0)
    {
      //std::cout << "step = " << ii << std::endl;
    }
    fun->setParam(x);
    double val = fun->f();
    Eigen::VectorXd grad = fun->df();
    double norm = infNorm(grad);
    if (ii==0)
    {
      h = maxStep / norm;
    }
    double init_h = h;
    int ret = lineSearch(fun, x, grad, h);
    std::cout << "step = " << ii << "                              " << " f = " << val << " norm = " << norm << " " << "h = " << init_h << " " << h << std::endl;
    if (ret < 0){
      break;
    }
    x0 = x;
    //oParameters.push_back(x0);
    oParameters.push_back(fun->distribution);
  }
}

int MaterialOptimizer::lineSearch(RealFun * fun, Eigen::VectorXd & x0, const Eigen::VectorXd & dir, double & h )
{
  double f0 = fun->f();
  Eigen::VectorXd grad0 = fun->df();
  double norm0 = infNorm(grad0);
  double norm1 = norm0;
  double f1 = f0;
  while (h>1.e-8){
    //minus sign here if we want to minimize function value.
    Eigen::VectorXd x = x0 - h*dir;
    clampVector(x, fun->lowerBounds, fun->upperBounds);
    fun->setParam(x);
    f1 = fun->f();
    Eigen::VectorXd grad1 = fun->df();
    norm1 = infNorm(grad1);
    //if gradient norm or function value decrease, return.
    if (/*norm1 < norm0 ||*/ f1 < f0){
      x0 = x;
      h *= 2;
      break;
    }
    else{
      h /= 2;
    }
  }
  if (fabs(norm0-norm1)<1.e-8 && fabs(f0-f1)<1.e-8)
  {
    std::cout << "diff_norm = " << fabs(norm0-norm1) << " diff_val = " << fabs(f0-f1) << std::endl;
    return -1;
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

Eigen::MatrixXd MaterialOptimizer::computeTargetStrains(int idim, float iYoungModulus, float iPoissonRatio, float iShearModulus)
{
  double E_x = iYoungModulus;
  double nu_xy = iPoissonRatio;
  double G_xy = iShearModulus;

  Eigen::MatrixXd G;
  // Cubic
  if (idim==2)
  {
    G = Eigen::MatrixXd::Identity(3,3);
    G(0,1) = -nu_xy;
    G(1,0) = -nu_xy;
    G *=  1./E_x;
    G(2,2) = 1./G_xy;
  }
  else
  {
    G = Eigen::MatrixXd::Identity(6,6);
    G(0,1) = -nu_xy;
    G(0,2) = -nu_xy;
    G(1,0) = -nu_xy;
    G(1,2) = -nu_xy;
    G(2,0) = -nu_xy;
    G(2,1) = -nu_xy;
    G *=  1./E_x;
    G(3,3) = 1./G_xy;
    G(4,4) = 1./G_xy;
    G(5,5) = 1./G_xy;
  }
  return G;
}

Eigen::MatrixXd MaterialOptimizer::computeTargetStrains(float iE_x, float iE_y, float iNu_xy, float iMu_xy)
{
  double E_x = iE_x;
  double E_y = iE_y;
  double nu_xy = iNu_xy;
  double G_xy = iMu_xy;

  // Cubic
  Eigen::MatrixXd G;
  G = Eigen::MatrixXd::Zero(3,3);
  G(0,0) = 1./E_x;
  G(0,1) = -nu_xy;
  G(1,0) = G(0,1);
  G(1,1) = 1./E_y;
  G(2,2) = 1./G_xy;
  return G;
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


