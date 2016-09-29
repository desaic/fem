#ifndef ExplorerUtilities_h
#define ExplorerUtilities_h

#include <vector>
#include <vtkVector.h>
#include "cfgDefs.h"

namespace ExplorerUtilities
{
  bool saveMicrostructure(const std::string &iFileName, int nx, int ny, const std::vector<int> &iMaterialAssignments);

  void getColorsForDistanceVisualization(const std::vector<cfgScalar> &iDistances, std::vector<vtkVector3i> &oColors);
  void getColorsForReducedCoordinatesVisualization(const std::vector<cfgScalar> &iPoints, int iDim, std::vector<vtkVector3i> &oColors);

  void writeMicrostructures(const std::string &iOutputDirectory, int iSize, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<cfgScalar> &iParameters, const std::string iSuffix="");
}

#endif

