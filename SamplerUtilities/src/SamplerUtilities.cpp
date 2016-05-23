#include "SamplerUtilities.h"

#include "qhullUtilities.h"

#include <Qhull.h>
#include <QhullFacetList.h>
#include <QhullVertex.h>
#include <QhullVertexSet.h>
using namespace orgQhull;

#include "DistanceTool.h"

#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

#include "cfgUtilities.h"
using namespace cfgUtil; 

void SamplerUtilities::getKMeans(int iNbIterations, int iNbClusters, const std::vector<cfgScalar> &iPoints, int iDim, std::vector<std::vector<int> > &oClusters, std::vector<int> *oCenters)
{
  oClusters.clear();

  std::vector<int> centers;
  int nbClusters = std::min(iNbClusters, (int)iPoints.size()/iDim);

  getFurthestPointsGreedy(nbClusters, iPoints, iDim, centers);
  std::vector<double> centerPositions = convertVec<cfgScalar,double>(getSubVector(iPoints, iDim, centers));

  std::vector<std::vector<int> > clusters;
  int iter;
  for (iter=0; iter<iNbIterations; iter++)
  {
    clusters.clear();
    clusters.resize(nbClusters);

    DistanceTool distanceTool(centerPositions);

    int ipoint, npoint=(int)iPoints.size()/iDim;
    for (ipoint=0; ipoint<npoint; ipoint++)
    {
      double p[3] = {iPoints[3*ipoint], iPoints[3*ipoint+1], iPoints[3*ipoint+2]};
      int closestCenter = distanceTool.getClosestPointIndex(p);
      clusters[closestCenter].push_back(ipoint);
    }
    centerPositions.clear();
    int icluster;
    for (icluster=0; icluster<nbClusters; icluster++)
    {
      std::vector<double> clusterPoints = convertVec<cfgScalar,double>(getSubVector(iPoints, iDim, clusters[icluster]));
      double center[3];
      cfgUtil::getBarycenter<double, 3>(clusterPoints, center);
      int icoord;
      for (icoord=0; icoord<3; icoord++)
      {
        centerPositions.push_back(center[icoord]);
      }
    }
  }
  int icluster;
  for (icluster=0; icluster<nbClusters; icluster++)
  {
    if (clusters[icluster].size()>0)
    {
      oClusters.push_back(clusters[icluster]);
    }
  }

  if (oCenters)
  {
    oCenters->clear();
    std::vector<double> points = convertVec<cfgScalar,double>(iPoints);
    DistanceTool distanceTool(points);
    int icluster;
    for (icluster=0; icluster<nbClusters; icluster++)
    {
      if (clusters[icluster].size()>0)
      {
        double p[3] = {centerPositions[3*icluster], centerPositions[3*icluster+1], centerPositions[3*icluster+2]};
        int closestPoint = distanceTool.getClosestPointIndex(p);
        oCenters->push_back(closestPoint);
      }
    }
  }
}

void SamplerUtilities::computeConvexHull(const std::vector<float> &iPoints, int iDim, std::vector<int> &oConvexHullVertices)
{
  assert(iPoints.size()%iDim==0);

  oConvexHullVertices.clear();
  coordT * points = cfgUtil::createArray<float,coordT>(iPoints);
  int numpoints = (int)iPoints.size()/iDim;
  boolT newIsMalloc = 0;

  /*qh_init_A (stdin, stdout, stderr, 0, NULL);
  int exitcode = setjmp (qh errexit);
  if (!exitcode)
  {
    qh_init_B(points, numpoints, iDim, newIsMalloc);
    qh_qhull();
    qh_check_output();
    print_summary();
  }*/ 

  std::string rboxCommand2, qhullCommand2;
  Qhull qhull(rboxCommand2.c_str(), iDim, numpoints, points, qhullCommand2.c_str());

  QhullVertexList vertices = qhull.vertexList();
  QhullLinkedList<QhullVertex>::const_iterator v_it, v_end=vertices.end();
  for (v_it=vertices.begin(); v_it!=v_end; ++v_it)
  {
    QhullPoint p = (*v_it).point();
    oConvexHullVertices.push_back(p.id());

    /*const coordT *coord = p.coordinates();
    int icoord;
    for (icoord=0; icoord<iDim; icoord++)
    {
      oConvexHullVertices.push_back(coord[icoord]);
    }*/ 
  }

  /*qhull.outputQhull();
  if(qhull.hasQhullMessage()){
    std::cout << "\nResults of qhull\n" << qhull.qhullMessage();
    qhull.clearQhullMessage();
  }

  QhullFacetList facets= qhull.facetList();
  std::cout << "\nFacets created by Qhull::runQhull()\n" << facets;

  QhullVertexList vertices = qhull.vertexList();
  std::cout << "\nVertices created by Qhull::runQhull()\n" << vertices; */ 
}

void SamplerUtilities::computeDelaundayTriangulation(const std::vector<float> &iPoints, int iDim, std::vector<int> &oFaceVertices, std::vector<int> &oBoundaryVertices, std::vector<int> &oBoundaryFaces, std::vector<float> *oDistancesToBoundary)
{
  assert(iPoints.size()%iDim==0);

  oFaceVertices.clear();
  coordT * points = cfgUtil::createArray<float,coordT>(iPoints);
  int numpoints = (int)iPoints.size()/iDim;
  boolT newIsMalloc = 0;

  std::string rboxCommand2, qhullCommand2="d QJ";
  Qhull qhull(rboxCommand2.c_str(), iDim, numpoints, points, qhullCommand2.c_str());

  std::vector<int> allFaceVertices;
  qhullUtilities::getFacetsVertices(qhull.facetList(), allFaceVertices);
  //oFaceVertices = allFaceVertices;

  Eigen::Vector3f meanLengths = getMeanEdgeLengthPerAxis(iPoints, allFaceVertices, 4);
  Eigen::Vector3f medianLengths = getMedianEdgeLengthPerAxis(iPoints, allFaceVertices, 4);

  float meanEdgeLength = getMeanEdgeLength(iPoints, allFaceVertices, 4);
  float minEdgeLength = getMinEdgeLength(iPoints, allFaceVertices, 4);
  //qhullUtilities::getBoundaryVertices(qhull, qhull.facetList(), 1.65, *oBoundaryVertices);

  std::vector<QhullFacet> smallFacets, largeFacets;
  //float maxEdgeSize = 1.8; //3D
  float maxEdgeSize = 2.5; //2D 
  //qhullUtilities::sortFaces(qhull, qhull.facetList(), maxEdgeSize, smallFacets, largeFacets);

  maxEdgeSize = 3; //10
  std::vector<float> scale;
  scale.push_back(1.f/meanLengths[0]);
  scale.push_back(1.f/meanLengths[1]);
  scale.push_back(1.f/meanLengths[2]);
  qhullUtilities::sortFaces(qhull, qhull.facetList(), scale, maxEdgeSize, smallFacets, largeFacets);

  qhullUtilities::getBoundaryVertices(smallFacets, oBoundaryVertices, oBoundaryFaces);

  oFaceVertices.clear();
  qhullUtilities::getFacetsVertices(smallFacets, oFaceVertices);

  /*oFaceVertices.clear();
  qhullUtilities::filterFaces(qhull, qhull.facetList(), 1.5, oFaceVertices);
  *oBoundaryVertices = oFaceVertices;*/ 

  if (oDistancesToBoundary)
  {
    oDistancesToBoundary->clear();

    QhullPoints points = qhull.points();
    int ipoint=0, npoint=(int)points.size();
    for (ipoint=0; ipoint<npoint; ipoint++)
    {
      QhullPoint p = points.at(ipoint);
      float dist = qhullUtilities::distToClosestFacet(p, smallFacets);
      oDistancesToBoundary->push_back(dist);
    }
  }

  if (0)
  {
    std::ofstream stream("Summary.txt");
    qhull.outputQhull();
    if(qhull.hasQhullMessage()){
      stream << "\nResults of qhull\n" << qhull.qhullMessage();
      qhull.clearQhullMessage();
    }

    QhullFacetList facets= qhull.facetList();
    stream << "\nFacets created by Qhull::runQhull()\n" << facets;

    QhullVertexList vertices = qhull.vertexList();
    stream << "\nVertices created by Qhull::runQhull()\n" << vertices;
  }
}

void SamplerUtilities::computeDistancesToPointCloud(const std::vector<float> &iPoints, const std::vector<float> &iPointCloud, std::vector<float> &oDistances)
{
  oDistances.clear();

  std::vector<double> cloudPositions = convertVec<float,double>(iPointCloud);
  std::vector<double> pointsToTest = convertVec<float,double>(iPoints);

  DistanceTool distanceTool(cloudPositions);
  int ipoint=0, npoint=(int)pointsToTest.size()/3;
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    double ClosestPoint[3];
    distanceTool.getClosestPoint(&pointsToTest[3*ipoint], ClosestPoint);
    float Dist = (Eigen::Vector3f(ClosestPoint[0], ClosestPoint[1], ClosestPoint[2]) - Eigen::Vector3f(pointsToTest[3 * ipoint], pointsToTest[3 * ipoint + 1], pointsToTest[3 * ipoint + 2])).norm();
    oDistances.push_back(Dist);
  }
}

void SamplerUtilities::getClosestPoints(const std::vector<float> &iPoints, int iDim, std::vector<int> &iRefPointIndices, float iRange, std::vector<int> &oPoints)
{
  oPoints.clear();
  std::vector<float> pointCloud = getSubVector(iPoints, iDim, iRefPointIndices);
  std::vector<float> distances;
  computeDistancesToPointCloud(iPoints, pointCloud, distances);
  int ipoint=0, npoint=(int)distances.size();
  for (ipoint=0; ipoint<npoint; ipoint++)
  {
    if (distances[ipoint] < iRange)
    {
      oPoints.push_back(ipoint);
    }
  }
}


void SamplerUtilities::samplePointsGreedyV2(int iOutputNbPoints, const std::vector<cfgScalar> &iPoints, int iDim, cfgScalar iMinRadius, std::vector<int> &oPointIndices)
{
  oPointIndices.clear();

  cfgScalar minRadiusSquare = iMinRadius*iMinRadius;

  std::set<int> points;
  int ipoint, npoint=(int)iPoints.size()/iDim;
  for (ipoint=1; ipoint<npoint; ipoint++)
  {
    points.insert(ipoint);
  }
  while (points.size()>0)
  {
    int indPoint = *points.begin();
    oPointIndices.push_back(indPoint);

    std::set<int>::iterator it, it_end=points.end();
    std::vector<int> pointsToErase;
    cfgScalar minDist = FLT_MAX;
    for (it=points.begin(); it!=it_end; it++)
    {
      int indVertex = *it;
      cfgScalar dist =0;
      for (int icoord=0; icoord<iDim; icoord++)
      {
        cfgScalar a = iPoints[indPoint*iDim+icoord];
        cfgScalar b = iPoints[indVertex*iDim+icoord];
        dist += (a-b)*(a-b);
      }

      if (dist < minRadiusSquare)
      {
        pointsToErase.push_back(indVertex);
      }
      if (dist < minDist)
      {
        minDist = dist;
      }
    }
    points.erase(indPoint);
    std::cout << "nb erase points: " << pointsToErase.size() << " minDist = " << minDist << std::endl;
    for (int ipoint=0; ipoint<pointsToErase.size(); ipoint++)
    {
      points.erase(pointsToErase[ipoint]);
    }
    std::cout << "npoints: " << points.size() << std::endl;
  }


  for (ipoint=1; ipoint<iOutputNbPoints && points.size()>0; ipoint++)
  {
    //if (ipoint % 100 == 1)
    {
      std::cout << "indPoint = " << ipoint << std::endl;
    }
    cfgScalar maxDist = 0;
    int furthestPoint = *points.begin();
    std::set<int>::iterator it, it_end=points.end();
    for (it=points.begin(); it!=it_end; it++)
    {
      int indVertex = *it;

      cfgScalar currentMinDist = FLT_MAX;
      int jpoint, njpoint=(int)oPointIndices.size();
      for (jpoint=0; jpoint<njpoint; jpoint++)
      {
        int indPoint2 = oPointIndices[jpoint];
        cfgScalar dist=0;
        for (int icoord=0; icoord<iDim; icoord++)
        {
          cfgScalar a = iPoints[indVertex*iDim+icoord];
          cfgScalar b = iPoints[indPoint2*iDim+icoord];
          dist += (a-b)*(a-b);
        }
        if (dist < currentMinDist)
        {
           currentMinDist = dist;
        }
      }
      if (currentMinDist > maxDist)
      {
        maxDist = currentMinDist;
        furthestPoint = indVertex;
      }
    }
    oPointIndices.push_back(furthestPoint);
    points.erase(furthestPoint);

    std::vector<int> pointsToErase;
    cfgScalar minDist = FLT_MAX;
    for (it=points.begin(); it!=it_end; it++)
    {
      int indVertex = *it;
      cfgScalar dist =0;
      for (int icoord=0; icoord<iDim; icoord++)
      {
        cfgScalar a = iPoints[furthestPoint*iDim+icoord];
        cfgScalar b = iPoints[indVertex*iDim+icoord];
        dist += (a-b)*(a-b);
      }

      if (dist < minRadiusSquare)
      {
        pointsToErase.push_back(indVertex);
      }
      if (dist < minDist)
      {
        minDist = dist;
      }
    }
    std::cout << "nb erase points: " << pointsToErase.size() << " minDist = " << minDist << std::endl;
    for (int ipoint=0; ipoint<pointsToErase.size(); ipoint++)
    {
      points.erase(pointsToErase[ipoint]);
    }
    std::cout << "npoints: " << points.size() << std::endl;
  }
}

void SamplerUtilities::samplePointsGreedy(int iOutputNbPoints, const std::vector<cfgScalar> &iPoints, int iDim, cfgScalar iMinRadius, std::vector<int> &oPointIndices)
{
  oPointIndices.clear();

  cfgScalar minRadiusSquare = iMinRadius*iMinRadius;
  std::set<int> points;
  int ipoint, npoint=(int)iPoints.size()/iDim;
  for (ipoint=1; ipoint<npoint; ipoint++)
  {
    points.insert(ipoint);
  }
  int pointIndex = 0;
  oPointIndices.push_back(pointIndex);
 
  for (ipoint=1; ipoint<iOutputNbPoints && points.size()>0; ipoint++)
  {
    //if (ipoint % 100 == 1)
    {
      std::cout << "indPoint = " << ipoint << std::endl;
    }
    cfgScalar maxDist = 0;
    int furthestPoint = *points.begin();
    std::set<int>::iterator it, it_end=points.end();
    for (it=points.begin(); it!=it_end; it++)
    {
      int indVertex = *it;

      cfgScalar currentMinDist = FLT_MAX;
      int jpoint, njpoint=(int)oPointIndices.size();
      for (jpoint=0; jpoint<njpoint; jpoint++)
      {
        int indPoint2 = oPointIndices[jpoint];
        cfgScalar dist=0;
        for (int icoord=0; icoord<iDim; icoord++)
        {
          cfgScalar a = iPoints[indVertex*iDim+icoord];
          cfgScalar b = iPoints[indPoint2*iDim+icoord];
          dist += (a-b)*(a-b);
        }
        if (dist < currentMinDist)
        {
           currentMinDist = dist;
        }
      }
      if (currentMinDist > maxDist)
      {
        maxDist = currentMinDist;
        furthestPoint = indVertex;
      }
    }
    oPointIndices.push_back(furthestPoint);
    points.erase(furthestPoint);

    std::vector<int> pointsToErase;
    cfgScalar minDist = FLT_MAX;
    for (it=points.begin(); it!=it_end; it++)
    {
      int indVertex = *it;
      cfgScalar dist =0;
      for (int icoord=0; icoord<iDim; icoord++)
      {
        cfgScalar a = iPoints[furthestPoint*iDim+icoord];
        cfgScalar b = iPoints[indVertex*iDim+icoord];
        dist += (a-b)*(a-b);
      }

      if (dist < minRadiusSquare)
      {
        pointsToErase.push_back(indVertex);
      }
      if (dist < minDist)
      {
        minDist = dist;
      }
    }
    std::cout << "nb erase points: " << pointsToErase.size() << " minDist = " << minDist << std::endl;
    for (int ipoint=0; ipoint<pointsToErase.size(); ipoint++)
    {
      points.erase(pointsToErase[ipoint]);
    }
    std::cout << "npoints: " << points.size() << std::endl;
  }
}

void SamplerUtilities::samplePointsRandom(int iOutputNbPoints, const std::vector<cfgScalar> &iPoints, int iDim, cfgScalar iMinRadius, std::vector<int> &oPointIndices)
{
  oPointIndices.clear();

  int npoint=(int)iPoints.size()/iDim;
  if (iOutputNbPoints < npoint)
  {
    cfgScalar minRadiusSquare = iMinRadius*iMinRadius;
    std::set<int> points;
    for (int ipoint=0; ipoint<npoint; ipoint++)
    {
      points.insert(ipoint);
    }

    for (int ipoint=0; ipoint<npoint-iOutputNbPoints; ipoint++)
    {
      if (points.size() % 10000 == 1)
      {
        std::cout << "npoints: " << points.size() << std::endl;
      }
      int indPoint = rand()%points.size();

      std::set<int>::iterator it = points.begin();
      for (int i=0; i<indPoint; i++)
      {
        it++;
      }
      points.erase(it);
    }
    oPointIndices = toStdVector(points);
  }
  else
  {
    oPointIndices = genIncrementalSequence(0, npoint-1);
  }
}

