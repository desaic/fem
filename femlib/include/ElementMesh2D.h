#ifndef ELEMENTMESH2D_HPP
#define ELEMENTMESH2D_HPP
#include <vector>
#include "vecmath.h"
#include "MatrixX.hpp"
#include "Eigen/Sparse"

class Element2D;
class Material2D;

class ElementMesh2D{
public:
  ///@brief elements will be deleted by destructor. Do not put same pointer in different meshes.
  ///Make copies if needed.
  std::vector<Element2D*>e;

  std::vector<Vector2f>x;
  
  ///@brief vertex positions at rest pose.
  std::vector<Vector2f>X;
  
  std::vector<Material2D*>m;
  
  ///@brief material for each element
  std::vector<int>me;
  
  ///@brief external forces applied to each dof.
  std::vector<Vector2f>fe;
  
  std::vector<int> fixed;

  std::vector<std::vector<Vector2f> > * u;
  ElementMesh2D();
  
  ///@brief utility. call after initializing or changing X and e 
  ///X is copied to x;
  void initArrays();
  ///@brief add a material to the list of materials in this mesh
  void addMaterial(Material2D*_m);
  ///@brief for debug, check the size of members.
  int check();

  float getEnergy();
  float getEnergy(int eIdx);

  ///@brief computes internal forces only
  std::vector<Vector2f> getForce(int eIdx);
  std::vector<Vector2f> getForce();

  MatrixXf getStiffness(int eIdx);

  ///@param trig if true, return only the upper triangle of the symmetric matrix.
  void getStiffnessSparse(std::vector<float> &val, bool trig = false, bool constrained=false, bool iFixedRigid=false);

  ///@param I row offsets. I.size() = matrix size + 1. I[size()-1]=number of non-zeros.
  ///@param J column indices.
  void stiffnessPattern(std::vector<int> & I, std::vector<int> & J, bool trig = false, bool iFixedRigid=false);

  MatrixXf getStiffness();
  
  float eleSize();

  virtual ~ElementMesh2D();

private:
  void fixRigid(Eigen::SparseMatrix<float> & K,  ElementMesh2D * mesh);

};


void getEleX(int ii, const ElementMesh2D * m, std::vector<Vector2f> &x);

void setEleX(int ii, ElementMesh2D * m, const std::vector<Vector2f> &x);

#endif
