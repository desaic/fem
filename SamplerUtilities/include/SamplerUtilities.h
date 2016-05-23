#ifndef SamplerUtilities_h
#define SamplerUtilities_h
#include <Eigen/Dense>
#include <string>
#include <map>

#include "cfgDefs.h"

namespace SamplerUtilities
{
  // Point cloud utilities
  // ----------------------
  void computeConvexHull(const std::vector<float> &iPoints, int iDim, std::vector<int> &oConvexHullVertices);
  void computeDelaundayTriangulation(const std::vector<float> &iPoints, int iDim, std::vector<int> &oFaces, std::vector<int> &oBoundaryVertices, std::vector<int> &oBoundaryFaces, std::vector<float> *oDistancesToBoundary=NULL);
  void getKMeans(int iNbIterations, int iNbClusters, const std::vector<cfgScalar> &iPoints, int iDim, std::vector<std::vector<int> > &oClusters, std::vector<int> *oCenters=NULL);
  
  void getClosestPoints(const std::vector<float> &iPoints, int iDim, std::vector<int> &iRefPointIndices, float iRange, std::vector<int> &oPoints);
  void computeDistancesToPointCloud(const std::vector<float> &iPoints, const std::vector<float> &iPointCloud, std::vector<float> &oDistances);

  void samplePointsGreedy(int iOutputNbPoints, const std::vector<cfgScalar> &iPoints, int iDim, cfgScalar iMinRadius, std::vector<int> &oPointIndices);
  void samplePointsGreedyV2(int iOutputNbPoints, const std::vector<cfgScalar> &iPoints, int iDim, cfgScalar iMinRadius, std::vector<int> &oPointIndices);
  void samplePointsRandom(int iOutputNbPoints, const std::vector<cfgScalar> &iPoints, int iDim, cfgScalar iMinRadius, std::vector<int> &oPointIndices);
};

#endif 

