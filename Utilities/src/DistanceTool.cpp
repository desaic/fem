/*=====================================================================================*/
/*! \file		DistanceTool.cpp
	\author		skourasm
	\brief		Implementation of class DistanceTool
 */
/*=====================================================================================*/

#include "DistanceTool.h"

void getPoint(const std::vector<double> &iPoints, int iPointIndex, int idim, double *ioPoint)
{
  for (int icoord=0; icoord<idim; icoord++)
  {
    ioPoint[icoord] = iPoints[idim*iPointIndex + icoord];
  }
}

double getSqrNorm(const double *iPoint1, const double *iPoint2, int idim)
{
  double SquareDist = 0;
  for (int icoord=0; icoord<idim; icoord++)
  {
    double Diff = iPoint1[icoord]-iPoint2[icoord];
    SquareDist += Diff*Diff;
  }
  return SquareDist;
}


DistanceTool::DistanceTool(int idim)
{
  assert(idim==2||idim==3);
  m_dim = idim;

  m_KDTree = NULL;
  m_positions = NULL;
}

DistanceTool::DistanceTool(const std::vector<double> &iX, int idim)
{
  assert(idim==2||idim==3);
  m_dim = idim;

  m_KDTree = NULL;
  m_positions = NULL;

  assert(m_X.size()%m_dim == 0);
  m_X = iX;

  int ivertex=0, nvertex=(int)iX.size()/m_dim;
  m_positions = new double[m_dim*nvertex];
  for (ivertex=0; ivertex<nvertex; ivertex++)
  {
    for (int icoord=0; icoord<m_dim; icoord++)
    {
      int ind = m_dim*ivertex+icoord;
      m_positions[ind] = iX[ind];
    }
  }
  m_KDTree = new KdTree(m_positions, m_dim, nvertex, 20);
}

DistanceTool::DistanceTool(const std::vector<double> &iX, const std::vector<int> &iIndices, int idim)
{
  assert(idim==2||idim==3);
  m_dim = idim;

  assert(m_X.size()%m_dim == 0);

  m_KDTree = NULL;
  m_positions = NULL;

  int ivertex=0, nvertex=(int)iIndices.size();
  for (ivertex=0; ivertex<nvertex; ivertex++)
  {
    for (int icoord=0; icoord<m_dim; icoord++)
    {
      int ind = m_dim*iIndices[ivertex]+icoord;
      m_X.push_back(iX[ind]);
    }
  }
  m_positions = new double[m_dim*nvertex];
  for (ivertex=0; ivertex<nvertex; ivertex++)
  {
    for (int icoord=0; icoord<m_dim; icoord++)
    {
      int ind = m_dim*iIndices[ivertex]+icoord;
      m_positions[m_dim*ivertex+icoord] = iX[ind];
    }
  }
  m_KDTree = new KdTree(m_positions, m_dim, nvertex, 20);
}

DistanceTool::~DistanceTool()
{
  if (m_KDTree) delete m_KDTree; m_KDTree = NULL;
  if (m_positions) delete[] m_positions; m_positions = NULL;
}

void DistanceTool::setPositions(const std::vector<double> &iX)
{
  assert(m_X.size()%m_dim == 0);
  m_X = iX;

  //m_KDTree.Init(iX);

  if (m_positions) delete[] m_positions;
  int ivertex=0, nvertex=(int)(iX.size()/m_dim);
  m_positions = new double[m_dim*nvertex];
  for (ivertex=0; ivertex<nvertex; ivertex++)
  {
    for (int icoord=0; icoord<m_dim; icoord++)
    {
      int ind = m_dim*ivertex+icoord;
      m_positions[ind] = iX[ind];
    }
  }
  if (m_KDTree) delete m_KDTree;
  m_KDTree = new KdTree(m_positions, m_dim, nvertex, 20);
}

int DistanceTool::getClosestPointIndex(const double *iPoint)
{
  int ClosestPointIndex = 0;

  if (0)
  {
    int ivertex=0, nvertex=(int)m_X.size()/m_dim;
    double MinSquareDist = DBL_MAX;
    for (ivertex=0; ivertex<nvertex; ivertex++)
    {
      double CurrentPoint[3];
      getPoint(m_X, ivertex, m_dim, CurrentPoint);
      double SquareDist = getSqrNorm(CurrentPoint, iPoint, m_dim);
      if (SquareDist < MinSquareDist)
      {
        MinSquareDist = SquareDist;
        ClosestPointIndex = ivertex;
      }
    }
  }
  else
  {
    m_KDTree->setNOfNeighbours(1);
    m_KDTree->queryPosition(iPoint);
    ClosestPointIndex = m_KDTree->getNeighbourPositionIndex(0);
  }
  return ClosestPointIndex;
}

void DistanceTool::getClosestPoint(const double *iPoint, double *ioClosestPoint)
{
  int ClosestPointIndex = getClosestPointIndex(iPoint);
  getPoint(m_X, ClosestPointIndex, m_dim, ioClosestPoint);
}

int DistanceTool::getClosestPointIndex(const double *iPoint, const std::set<int> &iPointsToIgnore)
{
  int ClosestPointIndex = 0;

  int ivertex=0, nvertex=(int)m_X.size()/m_dim;
  double MinSquareDist = DBL_MAX;
  for (ivertex=0; ivertex<nvertex; ivertex++)
  {
    if (iPointsToIgnore.count(ivertex)==0)
    {
      double CurrentPoint[3];
      getPoint(m_X, ivertex, m_dim, CurrentPoint);
      double SquareDist = getSqrNorm(CurrentPoint, iPoint, m_dim);
      if (SquareDist < MinSquareDist)
      {
        MinSquareDist = SquareDist;
        ClosestPointIndex = ivertex;
      }
    }
  }
  return ClosestPointIndex;
}

void DistanceTool::getClosestPoint(const double *iPoint, const std::set<int> &iPointsToIgnore, double *ioClosestPoint)
{
  int ClosestPointIndex = getClosestPointIndex(iPoint, iPointsToIgnore);
  getPoint(m_X, ClosestPointIndex, m_dim, ioClosestPoint);
}

void DistanceTool::getClosestPointIndices(const double *iPoint, double iRadius, const std::set<int> &iPointsToIgnore,
                              std::vector<int> &oClosestPoints)
{
  oClosestPoints.clear();

  if (0)
  {
    double SquareRadius = iRadius*iRadius;
    int ivertex=0, nvertex=(int)m_X.size()/m_dim;
    for (ivertex=0; ivertex<nvertex; ivertex++)
    {
      if (iPointsToIgnore.count(ivertex)==0)
      {
        double CurrentPoint[3];
        getPoint(m_X, ivertex, m_dim, CurrentPoint);
        double SquareDist = getSqrNorm(CurrentPoint, iPoint, m_dim);
        if (SquareDist < SquareRadius)
        {
          oClosestPoints.push_back(ivertex);
        }
      }
    }
  }
  else
  {
    m_KDTree->setNOfNeighbours((int)m_X.size()/m_dim);
    m_KDTree->queryRange(iPoint, iRadius*iRadius);
    int ipoint=0, npoint=m_KDTree->getNOfFoundNeighbours();
    for (ipoint=0; ipoint<npoint; ipoint++)
    {
      int indVertex = m_KDTree->getNeighbourPositionIndex(ipoint);
      if (iPointsToIgnore.count(indVertex)==0)
      {
        oClosestPoints.push_back(indVertex);
      }
    }
  }
}


