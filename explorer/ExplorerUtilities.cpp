#include "ExplorerUtilities.h"

#include <assert.h>
#include <fstream>
#include <iostream>

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



