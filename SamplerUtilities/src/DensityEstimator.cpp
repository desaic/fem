
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

cfgScalar DensityEstimator::evalKernel(const cfgScalar *x, const cfgScalar *xi)
{
  cfgScalar y = 0;
  for (int icoord=0; icoord<m_dim; icoord++)
  {
    cfgScalar diff = x[icoord]-xi[icoord];
    y += diff*diff;
  }
  y /= (m_radius*m_radius);
  cfgScalar t = 1-y;
  cfgScalar phi = t*t*t*t;
  return phi;
}

void DensityEstimator::getNeighbours(const cfgScalar *iP, std::vector<int> &oInds)
{
  assert(m_distanceTool);
  const std::set<int> pointsToIgnore;
  double * p = new double[m_dim];
  for (int icoord=0; icoord<m_dim; icoord++)
    p[icoord] = iP[icoord];
  m_distanceTool->getClosestPointIndices(&p[0], m_radius, pointsToIgnore, oInds);
  delete [] p;
}

cfgScalar DensityEstimator::computeDensity(const cfgScalar *iP)
{
  double density = 0;
  std::vector<int> neighbours;
  getNeighbours(iP, neighbours);
  int ivertex, nvertex=(int)neighbours.size();
  for (ivertex=0; ivertex<nvertex; ivertex++)
  {
    int indVertex = neighbours[ivertex];
    cfgScalar * Pi = &m_x[indVertex*m_dim];
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
    if (ipoint % 1000 == 0)
      std::cout << ipoint << "/" << npoint << std::endl;
    cfgScalar * P = &m_x[m_dim*ipoint];
    cfgScalar density = computeDensity(P);
    densities.push_back(density);
  }
  return densities;
}

std::vector<cfgScalar> DensityEstimator::computeDensities(int iIndexStart, int iIndexEnd, std::vector<cfgScalar> &ioPreviousDensities)
{
  std::set<int> touchedPoints;
  std::vector<cfgScalar> densities;
  for (int ipoint=iIndexStart; ipoint<iIndexEnd; ipoint++)
  {
    if (ipoint % 1000 == 0)
      std::cout << ipoint << "/" << iIndexEnd-iIndexStart << std::endl;

    cfgScalar * P = &m_x[m_dim*ipoint];
    cfgScalar density = computeDensity(P);
    densities.push_back(density);

    std::vector<int> neighbours;
    getNeighbours(P, neighbours);
    int nneighbour = (int)neighbours.size();
    for (int ineighbour=0; ineighbour<nneighbour; ineighbour++)
    {
      int indPoint = neighbours[ineighbour];
      if (indPoint<iIndexStart)
      {
        cfgScalar * Q = &m_x[m_dim*indPoint];
        double phi = evalKernel(Q, P);
        ioPreviousDensities[indPoint] += phi;
      }
    }
  }
  return densities;
}

std::vector<cfgScalar> DensityEstimator::computeDensities(const std::vector<int> &iPointIndices)
{
  std::vector<cfgScalar> densities;
  int npoint=(int)iPointIndices.size();
  for (int ipoint=0; ipoint<npoint; ipoint++)
  {
    if (ipoint % 1000 == 0)
      std::cout << ipoint << "/" << npoint << std::endl;

    int indPoint = iPointIndices[ipoint];
    cfgScalar * P = &m_x[m_dim*indPoint];
    cfgScalar density = computeDensity(P);
    densities.push_back(density);
  }
  return densities;
}






