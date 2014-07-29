#ifndef ELEMENTMESH_HPP
#define ELEMENTMESH_HPP
#include <vector>
#include "vecmath.h"
#include "BoundVec3.hpp"

class Element;
class Material;
class ElementMesh{
public:
  std::vector<Element*>e;
  std::vector<Vector3f>x;
  std::vector<Material*>m;
  ///@brief material for each element
  std::vector<int>me;
  ///@brief external forces applied to each dof.
  std::vector<Vector3f>fe;
  ///@brief vertex bounds.
  std::vector<BoundVec3>bd;
  ElementMesh();
  
  ///@brief for debug, check the size of members.
  int check();

  float getEnergy();
  float getEnergy(int eIdx);
};
#endif