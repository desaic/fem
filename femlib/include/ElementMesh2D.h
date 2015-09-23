#ifndef ELEMENTMESH2D_HPP
#define ELEMENTMESH2D_HPP
#include <vector>
#include "vecmath.h"
#include "MatrixX.hpp"
#include "Eigen/Sparse"
#include "cfgDefs.h"

class Element2D;
class Material2D;

class ElementMesh2D{
public:
  ///@brief elements will be deleted by destructor. Do not put same pointer in different meshes.
  ///Make copies if needed.
  std::vector<Element2D*>e;

  std::vector<Vector2S>x;
  
  ///@brief vertex positions at rest pose.
  std::vector<Vector2S>X;
  
  std::vector<Material2D*>m;
  
  ///@brief material for each element
  std::vector<int>me;
  
  ///@brief external forces applied to each dof.
  std::vector<Vector2S>fe;
  
  std::vector<int> fixed;

  std::vector<std::vector<Vector2S> > * u;
  ElementMesh2D();
  
  ///@brief utility. call after initializing or changing X and e 
  ///X is copied to x;
  void initArrays();
  ///@brief add a material to the list of materials in this mesh
  void addMaterial(Material2D*_m);
  ///@brief for debug, check the size of members.
  int check();

  cfgScalar getEnergy();
  cfgScalar getEnergy(int eIdx);

  ///@brief computes internal forces only
  std::vector<Vector2S> getForce(int eIdx);
  std::vector<Vector2S> getForce();

  MatrixXS getStiffness(int eIdx);

  ///@param trig if true, return only the upper triangle of the symmetric matrix.
  void getStiffnessSparse(std::vector<cfgScalar> &val, bool trig = false, bool constrained=false, bool iFixedTranslation=false, bool iFixedRotation=false, bool iPeriodic=false);


  ///@param I row offsets. I.size() = matrix size + 1. I[size()-1]=number of non-zeros.
  ///@param J column indices.
  void stiffnessPattern(std::vector<int> & I, std::vector<int> & J, bool trig = false, bool iFixedTranslation=false, bool iFixedRotation=false, bool iPeriodic=false);

  MatrixXS getStiffness();
  
  cfgScalar eleSize();

  virtual ~ElementMesh2D();

private:
  void fixRigid(Eigen::SparseMatrix<cfgScalar> & K, bool iTriangular, ElementMesh2D * mesh);
  void fixTranslation(Eigen::SparseMatrix<cfgScalar> & K, bool iTriangular, ElementMesh2D * mesh);
  void fixRotation(Eigen::SparseMatrix<cfgScalar> & K, bool iTriangular, ElementMesh2D * mesh);
  void enforcePeriodicity(Eigen::SparseMatrix<cfgScalar> & K, bool iTriangular, ElementMesh2D * mesh);
};


void getEleX(int ii, const ElementMesh2D * m, std::vector<Vector2S> &x);

void setEleX(int ii, ElementMesh2D * m, const std::vector<Vector2S> &x);

#endif
