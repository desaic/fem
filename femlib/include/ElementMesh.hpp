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
  ///@brief elements will be deleted by destructor. Do not put same pointer in different meshes.
  ///Make copies if needed.
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

  std::vector<std::vector<Vector3f> > * u;
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

  ///@param trig if true, return only the upper triangle of the symmetric matrix.
  void getStiffnessSparse(std::vector<float> &val, bool trig = false, bool constrained=false);

  ///@param I row offsets. I.size() = matrix size + 1. I[size()-1]=number of non-zeros.
  ///@param J column indices.
  void stiffnessPattern(std::vector<int> & I, std::vector<int> & J, bool trig = false);

  MatrixXd getStiffness();

  virtual ~ElementMesh();
};


void getEleX(int ii, const ElementMesh * m, std::vector<Vector3f> &x);

void setEleX(int ii, ElementMesh * m, const std::vector<Vector3f> &x);

#endif
