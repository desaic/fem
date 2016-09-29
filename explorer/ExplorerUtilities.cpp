#include "ExplorerUtilities.h"

#include <assert.h>
#include <fstream>
#include <iostream>
#include "cfgDefs.h"
#include <cfgMaterialUtilities.h>
#include <cfgUtilities.h>

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

void ExplorerUtilities::getColorsForReducedCoordinatesVisualization(const std::vector<cfgScalar> &iPoints, int iDim, std::vector<vtkVector3i> &oColors)
{
  oColors.clear();

  assert(iDim==3);
  std::vector<cfgScalar> points = iPoints;

  std::vector<cfgScalar> lengths(iDim, 1);
  bool resOk = cfgMaterialUtilities::rescaleData(points, iDim, lengths);
  int npoint = (int)points.size()/iDim;
  for (int ipoint=0; ipoint<npoint; ipoint++)
  {
    vtkVector3i color(points[3*ipoint]*255, points[3*ipoint+1]*255, points[3*ipoint+2]*255);
    oColors.push_back(color);
  }
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

  vtkVector3i colMin(0, 0, 255);
  vtkVector3i colMax(255, 0, 0);

  vtkVector3i black(0, 0, 0);
  vtkVector3i colZero(255,255,255);

  colMin = black;
  colMax = vtkVector3i(220, 220, 220);
  colZero = colMax;

  cfgScalar coeff = 0.1;
  std::vector<vtkVector3i> colors;
  for (int ipoint=0; ipoint<npoint; ipoint++)
  {
    vtkVector3i matColor;
    if (iDistances[ipoint] < 0)
    {
      cfgScalar w = iDistances[ipoint]/minValue;
      matColor = vtkVector3i( w*colMin[0] + (1-w)*colZero[0], w*colMin[1] + (1-w)*colZero[1], w*colMin[2] + (1-w)*colZero[2]);
    }
    else if (iDistances[ipoint] > 0)
    {
      cfgScalar w = iDistances[ipoint]/maxValue;
      matColor = vtkVector3i(w*255 + (1-w)*255, (1-w)*255, (1-w)*255);
      matColor = vtkVector3i( w*colMax[0] + (1-w)*colZero[0], w*colMax[1] + (1-w)*colZero[1], w*colMax[2] + (1-w)*colZero[2]);
    }
    else
    {
      matColor = vtkVector3i(255, 0, 0);
    }
    matColor.Set((1-coeff)*matColor[0]+coeff*black[0], (1-coeff)*matColor[1]+coeff*black[1], (1-coeff)*matColor[2]+coeff*black[2]);
    oColors.push_back(matColor);
  }
}

void ExplorerUtilities::writeMicrostructures(const std::string &iOutputDirectory, int iSize, const std::vector<std::vector<int> > &iMaterialAssignments, const std::vector<cfgScalar> &iParameters, const std::string iSuffix)
{
  std::string fileRootName = iOutputDirectory + "level" + std::to_string(iSize) + "_";
  if (iSuffix!="")
  {
    fileRootName += iSuffix + "_";
  }
  std::string fileExtension = ".bin";
  cfgUtil::writeBinary<float>(fileRootName + "params" + fileExtension, iParameters);
  cfgUtil::writeBinary<int>(fileRootName + "matAssignments" + fileExtension, iMaterialAssignments);
}

