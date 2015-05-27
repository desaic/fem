#include "qHullUtilities.h"

#include <QhullVertex.h>
#include <QhullVertexSet.h>
#include <QhullSets.h>
#include <QhullRidge.h>
#include <QhullFacetSet.h>
#include <QhullPoint.h>
#include <Qhull.h>

#include <Vector3f.h>
#include <set>

#include "cfgUtilities.h"
#include "cfgMaterialUtilities.h"
using namespace cfgMaterialUtilities;

void qhullUtilities::filterFaces(const Qhull &iHull, const QhullFacetList &iFacets, float iMaxEdgeLength, std::vector<int> &oFaceVertices)
{
  std::vector<QhullFacet> smallFacets, largeFacets;
  sortFaces(iHull, iFacets, iMaxEdgeLength, smallFacets, largeFacets);
  getFacetsVertices(smallFacets, oFaceVertices);
  return ;

  QhullPoints points = iHull.points();

  QhullFacetList::const_iterator f_it, f_end=iFacets.end();
  for (f_it=iFacets.begin(); f_it!=f_end; ++f_it)
  {
    if ((*f_it).isGood())
    {
      std::vector<std::vector<int> > ridges;
      getRidges((*f_it), ridges);

      bool longEdgeFound = false;

      int iridge=0, nridge=(int)ridges.size();
      for (iridge=0; iridge<nridge; iridge++)
      {
        const std::vector<int> & ridgeVertices = ridges[iridge];
        int ivertex=0, nvertex=(int)ridgeVertices.size();
        for (ivertex=0; ivertex<nvertex; ivertex++)
        {
          int indP1 = ridgeVertices[ivertex];
          int indP2 = ridgeVertices[(ivertex+1)%nvertex];

          QhullPoint p1 = points.at(indP1);
          QhullPoint p2 = points.at(indP2);

          Vector3f q1(p1[0],p1[1],p1[2]);
          Vector3f q2(p2[0],p2[1],p2[2]);

          float length = (q1-q2).abs();
          if (length>iMaxEdgeLength)
          {
            longEdgeFound = true;

            /*oFaceVertices.push_back(indP1);
            oFaceVertices.push_back(indP2);

            oFaceVertices.push_back(indP1);
            oFaceVertices.push_back(indP2);*/ 
          }
        }
      }
      if (longEdgeFound)
      {
        addFacetVertices((*f_it),oFaceVertices);
      }
    }
  }
}

void qhullUtilities::sortFaces(const Qhull &iHull, const QhullFacetList &iFacets, float iMaxEdgeLength, std::vector<QhullFacet> &oSmallFacets, std::vector<QhullFacet> &oLargeFacets)
{
  oSmallFacets.clear();
  oLargeFacets.clear();

  QhullPoints points = iHull.points();
  QhullFacetList::const_iterator f_it, f_end=iFacets.end();
  for (f_it=iFacets.begin(); f_it!=f_end; ++f_it)
  {
    if ((*f_it).isGood())
    {
      std::vector<std::vector<int> > ridges;
      getRidges((*f_it), ridges);

      bool longEdgeFound = false;

      int iridge=0, nridge=(int)ridges.size();
      for (iridge=0; iridge<nridge; iridge++)
      {
        const std::vector<int> & ridgeVertices = ridges[iridge];
        int ivertex=0, nvertex=(int)ridgeVertices.size();
        for (ivertex=0; ivertex<nvertex; ivertex++)
        {
          int indP1 = ridgeVertices[ivertex];
          int indP2 = ridgeVertices[(ivertex+1)%nvertex];

          QhullPoint p1 = points.at(indP1);
          QhullPoint p2 = points.at(indP2);

          Vector3f q1(p1[0],p1[1],p1[2]);
          Vector3f q2(p2[0],p2[1],p2[2]);

          float length = (q1-q2).abs();
          if (length>iMaxEdgeLength)
          {
            longEdgeFound = true;

            /*oFaceVertices.push_back(indP1);
            oFaceVertices.push_back(indP2);

            oFaceVertices.push_back(indP1);
            oFaceVertices.push_back(indP2);*/ 
          }
        }
      }
      if (longEdgeFound)
      {
        oLargeFacets.push_back(*f_it);
      }
      else
      {
        oSmallFacets.push_back(*f_it);
      }
    }
  }
}

void qhullUtilities::addFacetVertices(const QhullFacet &iFacet, std::vector<int> &ioFaceVertices)
{
  QhullVertexSet fv = iFacet.vertices();
  if (!fv.isEmpty())
  {
    QhullVertexSetIterator fv_it = fv;
    while(fv_it.hasNext())
    {
      QhullVertex v = fv_it.next();
      QhullPoint p = v.point();
      ioFaceVertices.push_back(p.id());      
    }
  }
}

void qhullUtilities::getFacetsVertices(const QhullFacetList &iFacets, std::vector<int> &oFaceVertices)
{
  oFaceVertices.clear();
  QhullFacetList::const_iterator f_it, f_end=iFacets.end();
  for (f_it=iFacets.begin(); f_it!=f_end; ++f_it)
  {
    if ((*f_it).isGood())
    {
      addFacetVertices(*f_it, oFaceVertices);
    }
  }
}

void qhullUtilities::getFacetsVertices(const std::vector<QhullFacet> &iFacets, std::vector<int> &oFaceVertices)
{
  int iface=0, nface=(int)iFacets.size();
  for (iface=0; iface<nface; iface++)
  {
    QhullFacet f = iFacets[iface];
    if (f.isGood())
    {
      addFacetVertices(f, oFaceVertices);
    }
  }
}

void qhullUtilities::getBoundaryVertices(const std::vector<QhullFacet> &iFacets, std::vector<int> &oBoundaryVertices, std::vector<int> &oBoundaryFaces)
{
  std::set<int> FacetSet;
  int ifacet=0, nfacet=(int)iFacets.size();
  for (ifacet=0; ifacet<nfacet; ifacet++)
  {
    FacetSet.insert(iFacets[ifacet].id());
  }

  std::vector<QhullFacet> boundaryFacets;
  std::vector<std::vector<QhullFacet> > adjacentFacets;
  std::vector<QhullFacet>::const_iterator f_it, f_end=iFacets.end();
  for (f_it=iFacets.begin(); f_it!=f_end; f_it++)
  {
    QhullFacetSet ff = (*f_it).neighborFacets();
    if (!ff.isEmpty())
    {
       std::vector<QhullFacet> adjFacets;

      QhullFacetSetIterator ff_it = ff;
      while(ff_it.hasNext())
      {
        QhullFacet f = ff_it.next();
        if (FacetSet.count(f.id())>0)
        {
          adjFacets.push_back(f);
        }
      }
      if (adjFacets.size()<4)
      {
        boundaryFacets.push_back((*f_it));
        adjacentFacets.push_back(adjFacets);
      }
    }
  }

  std::set<int> BoundaryVertices;
  std::vector<int> BoundaryFaces;
  int iface=0, nface=(int)boundaryFacets.size();
  for (iface=0; iface<nface; iface++)
  {
    getBoundaryVertices(boundaryFacets[iface],adjacentFacets[iface], BoundaryVertices, BoundaryFaces);
  }
  oBoundaryVertices = cfgUtil::toStdVector(BoundaryVertices);
  oBoundaryFaces = BoundaryFaces;
 // getFacetsVertices(boundaryFacets, oBoundaryVertices);
} 

void qhullUtilities::getRidges(const QhullFacet &iFacet, std::vector<std::vector<int> > &oRidges)
{
  int tetfaces[4][3] = { {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 3, 0} };

  oRidges.clear();
  oRidges.resize(4);

  std::vector<int> facetVertices;
  addFacetVertices(iFacet, facetVertices);

  int iface=0;
  for (iface=0; iface<4; iface++)
  {
    int ipoint=0;
    for (ipoint=0; ipoint<4; ipoint++)
    {
      int indVertex = tetfaces[iface][ipoint%3];
      oRidges[iface].push_back(facetVertices[indVertex]);
    }
  }
}

void qhullUtilities::getRidgeVertices(const QhullFacet &iFacet, std::vector<std::set<int> > &oRidgeVertices)
{
  int tetfaces[4][3] = { {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 3, 0} };

  oRidgeVertices.clear();
  oRidgeVertices.resize(4);

  std::vector<int> facetVertices;
  addFacetVertices(iFacet, facetVertices);

  int iface=0;
  for (iface=0; iface<4; iface++)
  {
    int ipoint=0;
    for (ipoint=0; ipoint<3; ipoint++)
    {
      int indVertex = tetfaces[iface][ipoint];
      oRidgeVertices[iface].insert(facetVertices[indVertex]);
    }
  }
}

void qhullUtilities::getRidgeOrderedVertices(const QhullFacet &iFacet, std::vector<std::vector<int> > &oRidgeVertices)
{
  int tetfaces[4][3] = { {0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {3, 2, 0} };

  oRidgeVertices.clear();
  oRidgeVertices.resize(4);

  std::vector<int> facetVertices;
  addFacetVertices(iFacet, facetVertices);

  int iface=0;
  for (iface=0; iface<4; iface++)
  {
    int ipoint=0;
    for (ipoint=0; ipoint<3; ipoint++)
    {
      int indVertex = tetfaces[iface][ipoint];
      oRidgeVertices[iface].push_back(facetVertices[indVertex]);
    }
  }
}

void qhullUtilities::getBoundaryVertices(const QhullFacet &iFacet, const std::vector<QhullFacet> &iAdjFacets, std::set<int> &ioBoundaryVertices, std::vector<int> &ioBoundaryFaces)
{
  std::vector<std::set<int> > RidgeVertices;
  getRidgeVertices(iFacet, RidgeVertices);

  //std::vector<std::vector<int> > RidgeOrderedVertices;
  //getRidgeOrderedVertices(iFacet, RidgeOrderedVertices);

  std::vector<int> vecfaces;
  std::vector<int> adjFaceIndices(iAdjFacets.size());
  int iface=0, nface=(int)iAdjFacets.size();
  for (iface=0; iface<nface; iface++)
  {
    std::vector<std::set<int> > AdjRidgeVertices;
    getRidgeVertices(iAdjFacets[iface], AdjRidgeVertices);

    int index = -1;
    int iridge=0, nridge=(int)RidgeVertices.size();
    for (iridge=0; iridge<nridge && index<0; iridge++)
    {
      int iadjridge=0, nadjridge=(int)AdjRidgeVertices.size();
      for (iadjridge=0; iadjridge<nadjridge && index<0; iadjridge++)
      {
        if (RidgeVertices[iridge]==AdjRidgeVertices[iadjridge])
        {
          index = iridge;
        }
      }
    }
    adjFaceIndices[iface] = index;
  }

  std::vector<bool> hasAdjFacets(RidgeVertices.size(), false);
  for (iface=0; iface<nface; iface++)
  {
    int indAdjFace = adjFaceIndices[iface];
    hasAdjFacets[indAdjFace] = true;
  }

  int iridge=0, nridge=(int)RidgeVertices.size();
  for (iridge=0; iridge<nridge; iridge++)
  {
    if (!hasAdjFacets[iridge])
    {
      std::set<int>::const_iterator r_it, r_end=RidgeVertices[iridge].end();
      for (r_it=RidgeVertices[iridge].begin(); r_it!=r_end; r_it++)
      {
        ioBoundaryVertices.insert(*r_it);
        ioBoundaryFaces.push_back(*r_it);
      }
      //ioBoundaryFaces.insert(ioBoundaryFaces.end(), RidgeOrderedVertices[iridge].begin(), RidgeOrderedVertices[iridge].end());
    }
  }
}



