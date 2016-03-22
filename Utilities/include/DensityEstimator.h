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

  cfgScalar computeDensity(const cfgScalar *iP);
  std::vector<cfgScalar> computeDensities();
  std::vector<cfgScalar> computeDensities(int iIndexStart, int iIndexEnd, std::vector<cfgScalar> &ioPreviousDensities);
  std::vector<cfgScalar> computeDensities(const std::vector<int> &iPointIndices);

  cfgScalar estimateRadius();

private:
  cfgScalar evalKernel(const cfgScalar *x, const cfgScalar *xi);
  void getNeighbours(const cfgScalar *iP, std::vector<int> &oInds);

private:
  int m_dim;
  cfgScalar m_radius;
  std::vector<cfgScalar> m_x;
  DistanceTool * m_distanceTool;
};

#endif





