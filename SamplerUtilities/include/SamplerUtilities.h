#ifndef SamplerUtilities_h
#define SamplerUtilities_h
#include <Eigen/Dense>
#include <string>
#include <map>

#include "cfgDefs.h"

namespace SamplerUtilities
{
  // Point cloud utitlities
  // ----------------------
  void computeConvexHull(const std::vector<float> &iPoints, int iDim, std::vector<int> &oConvexHullVertices);
  void computeDelaundayTriangulation(const std::vector<float> &iPoints, int iDim, std::vector<int> &oFaces, std::vector<int> &oBoundaryVertices, std::vector<int> &oBoundaryFaces, std::vector<float> *oDistancesToBoundary=NULL);
  void getKMeans(int iNbIterations, int iNbClusters, const std::vector<cfgScalar> &iPoints, int iDim, std::vector<std::vector<int> > &oClusters, std::vector<int> *oCenters=NULL);
  
  void getClosestPoints(const std::vector<float> &iPoints, int iDim, std::vector<int> &iRefPointIndices, float iRange, std::vector<int> &oPoints);
  void computeDistancesToPointCloud(const std::vector<float> &iPoints, const std::vector<float> &iPointCloud, std::vector<float> &oDistances);
};

#endif 

