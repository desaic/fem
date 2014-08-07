#ifndef ELEMENTMESH_HPP
#define ELEMENTMESH_HPP
#include <vector>
#include "vecmath.h"
#include "BoundVec3.hpp"

class Element;
class Material;
struct MatrixXd;

class ElementMesh{
public:
  std::vector<Element*>e;
  std::vector<Vector3f>x;
  ///@brief vertex positions at rest pose.
  std::vector<Vector3f>X;
  std::vector<Material*>m;
  ///@brief material for each element
  std::vector<int>me;
  ///@brief external forces applied to each dof.
  std::vector<Vector3f>fe;
  std::vector<int> fixed;
  ElementMesh();
  
  ///@brief utility. call after initializing or changing X and e 
  ///X is copied to x;
  void initArrays();
  ///@brief for debug, check the size of members.
  int check();

  float getEnergy();
  float getEnergy(int eIdx);

  ///@brief computes internal forces only
  std::vector<Vector3f> getForce(int eIdx);
  std::vector<Vector3f> getForce();

  MatrixXd getStiffness(int eIdx);
  void getStiffnessSparse(std::vector<int> & I, std::vector<int> & J, std::vector<float>
    &val);

  MatrixXd getStiffness();
};
#endif