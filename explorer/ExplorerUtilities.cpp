#include "ExplorerUtilities.h"

#include <assert.h>
#include <fstream>
#include <iostream>
#include "cfgDefs.h"

bool ExplorerUtilities::saveMicrostructure(const std::string &iFileName, int nx, int ny, const std::vector<int> &iMaterialAssignments)
{
  assert(iMaterialAssignments.size()==nx*ny);


  std::ofstream stream(iFileName);
  if (!stream.is_open())
  {
    return false;
  }

  stream << nx << " " << ny << std::endl;
  int indMat = 0;
  for (int i=0; i<nx; i++)
  {
    for (int j=0; j<ny; j++)
    {
      int mat = iMaterialAssignments[indMat++];
      stream << mat << " ";
    }
    stream << std::endl;
  }
  return true;
}

void ExplorerUtilities::getColorsForDistanceVisualization(const std::vector<cfgScalar> &iDistances, std::vector<vtkVector3i> &oColors)
{
  oColors.clear();

  int npoint = (int)iDistances.size();

  cfgScalar minValue = FLT_MAX;
  cfgScalar maxValue = -FLT_MAX;
  for (int ipoint=0; ipoint<npoint; ipoint++)
  {
    cfgScalar dist = iDistances[ipoint];
    if (dist < minValue)
      minValue = dist;
    if (dist > maxValue)
      maxValue = dist;
  }
  std::cout << "min dist = " << minValue << ", max dist = " << maxValue << std::endl;

  vtkVector3i black(0,0,0);
  cfgScalar coeff = 0.1;
  std::vector<vtkVector3i> colors;
  for (int ipoint=0; ipoint<npoint; ipoint++)
  {
    vtkVector3i matColor;
    if (iDistances[ipoint] < 0)
    {
      cfgScalar w = iDistances[ipoint]/minValue;
      matColor = vtkVector3i((1-w)*255, (1-w)*255, w*255 + (1-w)*255);
    }
    else if (iDistances[ipoint] > 0)
    {
      cfgScalar w = iDistances[ipoint]/maxValue;
      matColor = vtkVector3i(w*255 + (1-w)*255, (1-w)*255, (1-w)*255);
    }
    else
    {
      matColor = vtkVector3i(255, 0, 0);
    }
    matColor.Set((1-coeff)*matColor[0]+coeff*black[0], (1-coeff)*matColor[1]+coeff*black[1], (1-coeff)*matColor[2]+coeff*black[2]);
    oColors.push_back(matColor);
  }
}


