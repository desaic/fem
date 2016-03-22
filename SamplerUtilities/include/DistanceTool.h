/*=====================================================================================*/
/*! \file		DistanceTool.h
	\author		skourasm
	\brief		Declaration of class DistanceTool
 */
/*=====================================================================================*/

#ifndef DistanceTool_h
#define DistanceTool_h

#include "KdTree.h"
#include <set>

class DistanceTool 
{
public: 
  DistanceTool(int idim=3);
  DistanceTool(const std::vector<double> &iX, int idim=3);
  DistanceTool(const std::vector<double> &iX, const std::vector<int> &iIndices, int idim=3);

  virtual ~DistanceTool();

  void setPositions(const std::vector<double> &iX);
  const std::vector<double> & setPositions() { return m_X;};

  int getClosestPointIndex(const double *iPoint);
  void getClosestPoint(const double *iPoint, double *ioClosestPoint);

  int getClosestPointIndex(const double *iPoint, const std::set<int> &iPointsToIgnore);
  void getClosestPoint(const double *iPoint, const std::set<int> &iPointsToIgnore, double *ioClosestPoint);

  void getClosestPointIndices(const double *iPoint, double iRadius, const std::set<int> &iPointsToIgnore,
                              std::vector<int> &oClosestPoints);

protected:
  int m_dim;
  std::vector<double> m_X;

  KdTree *m_KDTree;
  double * m_positions;
};

#endif
