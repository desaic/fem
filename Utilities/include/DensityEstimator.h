#ifndef DensityEstimator_h
#define DensityEstimator_h

#include "cfgDefs.h"

#include "DistanceTool.h"

class DensityEstimator
{
public:
  DensityEstimator(int iDim);
  ~DensityEstimator();

  void setKernelRadius(cfgScalar iRadius); 
  void setPoints(const std::vector<cfgScalar> &iPoints);
  void init();

  cfgScalar computeDensity(const Vector3S &iP);
  std::vector<cfgScalar> computeDensities();

  cfgScalar estimateRadius();

private:
  cfgScalar evalKernel(const Vector3S &x, const Vector3S &xi);
  void getNeighbours(const Vector3S &iP, std::vector<int> &oInds);

private:
  int m_dim;
  cfgScalar m_radius;
  std::vector<cfgScalar> m_x;
  DistanceTool * m_distanceTool;
};

#endif





