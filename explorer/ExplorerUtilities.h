#ifndef ExplorerUtilities_h
#define ExplorerUtilities_h

#include <vector>

namespace ExplorerUtilities
{
  bool saveMicrostructure(const std::string &iFileName, int nx, int ny, const std::vector<int> &iMaterialAssignments);
}

#endif

