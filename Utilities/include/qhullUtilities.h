#ifndef qhullUtilities_h
#define qhullUtilities_h

#include <QhullFacet.h>
#include <QhullFacetList.h>
using namespace orgQhull;

#include <set>

namespace qhullUtilities
{
  void addFacetVertices(const QhullFacet &iFacet, std::vector<int> &ioFaceVertices);
  void getFacetsVertices(const QhullFacetList &iFacets, std::vector<int> &oFaceVertices);
  void getFacetsVertices(const std::vector<QhullFacet> &iFacets, std::vector<int> &oFaceVertices);

  void filterFaces(const Qhull &iHull, const QhullFacetList &iFacets, float iMaxEdgeLength, std::vector<int> &oFaceVertices);
  void sortFaces(const Qhull &iHull, const QhullFacetList &iFacets, float iMaxEdgeLength, std::vector<QhullFacet> &oSmallFacets, std::vector<QhullFacet> &oLargeFacets);
  void sortFaces(const Qhull &iHull, const QhullFacetList &iFacets, const std::vector<float> &iScale, float iMaxEdgeLength, std::vector<QhullFacet> &oSmallFacets, std::vector<QhullFacet> &oLargeFacets);

  void getBoundaryVertices(const QhullFacet &iFacet, const std::vector<QhullFacet> &iAdjFacets, std::set<int> &ioBoundaryVertices, std::vector<int> &ioBoundaryFaces);
  void getBoundaryVertices(const std::vector<QhullFacet> &iFacets, std::vector<int> &oBoundaryVertices, std::vector<int> &oBoundaryFaces);

  void getRidgeVertices(const QhullFacet &iFacet, std::vector<std::set<int> > &oRidgeVertices);
  void getRidgeOrderedVertices(const QhullFacet &iFacet, std::vector<std::vector<int> > &oRidgeVertices);
  void getRidges(const QhullFacet &iFacet, std::vector<std::vector<int> > &oRidges);

  float distToClosestFacet(QhullPoint &iPoint, const std::vector<QhullFacet> &iFacets);
};

#endif 