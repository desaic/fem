
#include "DensityEstimator.h"

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

#include "cfgUtilities.h."
using namespace cfgUtil;

DensityEstimator::DensityEstimator(int iDim)
{
  m_dim = iDim;
  m_radius = 0;
  m_distanceTool = NULL;
}

DensityEstimator::~DensityEstimator()
{
  SAFE_DELETE(m_distanceTool);
}

void DensityEstimator::setPoints(const std::vector<cfgScalar> &iPoints)
{
  m_x = iPoints;
}

void DensityEstimator::setKernelRadius(cfgScalar iRadius)
{
  m_radius = iRadius;
}

void DensityEstimator::init()
{
  if (m_radius==0)
  {
    m_radius = estimateRadius();
  }
  SAFE_DELETE(m_distanceTool);
  m_distanceTool = new DistanceTool(convertVec<cfgScalar, double>(m_x), m_dim);
}

cfgScalar DensityEstimator::estimateRadius()
{
  int nsubdiv = 10;
  std::vector<cfgScalar> box[2];
  getBoundingBox(m_x, m_dim, box);
  cfgScalar radius = FLT_MAX;
  for (int idim=0; idim<m_dim; idim++)
  {
    cfgScalar diff = box[1][idim]-box[0][idim];
    radius = std::min(radius, diff/nsubdiv);
  }
  return radius;
}

cfgScalar DensityEstimator::evalKernel(const Vector3S &x, const Vector3S &xi)
{
  cfgScalar y = (x-xi).squaredNorm()/(m_radius*m_radius);
  cfgScalar t = 1-y;
  cfgScalar phi = t*t*t*t;
  return phi;
}

void DensityEstimator::getNeighbours(const Vector3S &iP, std::vector<int> &oInds)
{
  assert(m_distanceTool);
  const std::set<int> pointsToIgnore;
  double p[3] = {iP[0],iP[1],iP[2]};
  m_distanceTool->getClosestPointIndices(&p[0], m_radius, pointsToIgnore, oInds);
}

cfgScalar DensityEstimator::computeDensity(const Vector3S &iP)
{
  double density = 0;
  std::vector<int> neighbours;
  getNeighbours(iP, neighbours);
  int ivertex, nvertex=(int)neighbours.size();
  for (ivertex=0; ivertex<nvertex; ivertex++)
  {
    int indVertex = neighbours[ivertex];
    Vector3S Pi = getVector3S(indVertex, m_x);
    double phi = evalKernel(iP, Pi);
    density += phi;
  }
  return (cfgScalar)density;
}

std::vector<cfgScalar> DensityEstimator::computeDensities()
{
  std::vector<cfgScalar> densities;
  int ipoint=0, npoint=(int)m_x.size()/m_dim;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    Vector3S P = getVector3S(ipoint, m_x);
    cfgScalar density = computeDensity(P);
    densities.push_back(density);
  }
  return densities;
}






