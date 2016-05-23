#ifndef ELEMENTMESH_HPP
#define ELEMENTMESH_HPP
#include <vector>
#include <Eigen/Sparse>
#include "cfgDefs.h"

class Element;
class Material;

class ElementMesh{
public:
  ///@brief elements will be deleted by destructor. Do not put same pointer in different meshes.
  ///Make copies if needed.
  std::vector<Element*>e;

  ///@brief deformed positions in t-1. used by dynamics(unit m)
//  std::vector<Vector3S>x0;

  ///@brief deformed positions (unit m)
  std::vector<Vector3S>x;
  
  ///@brief vertex positions at rest pose.
  std::vector<Vector3S>X;

  ///@brief velocity
  std::vector<Vector3S>v;

  ///@brief time step (unit s)
  cfgScalar dt;
  ///@brief gravitational constant
  Vector3S G;
  ///@TODO: change to per element density.
  ///@brief kg/m3
  cfgScalar density;
  ///@brief lumped mass
  std::vector<cfgScalar>mass;

  std::vector<Material*>m;
  
  ///@brief material for each element
  std::vector<int>me;
  
  ///@brief external forces applied to each dof.
  std::vector<Vector3S>fe;
  
  std::vector<int> fixed;
  
  ///@brief default constructor
  ElementMesh();
  
  ///@brief load a plain text hex fem mesh
  void load(std::istream & in, cfgScalar scale=1.0f);

  ///@brief utility. call after initializing or changing X and e 
  ///X is copied to x;
  void initArrays();

  ///@brief add a material to the list of materials in this mesh
  void addMaterial(Material*_m);

  ///@brief for debug, check the size of members.
  int check();

  cfgScalar getEnergy();
  cfgScalar getEnergy(int eIdx);

  ///@brief computes internal forces only
  std::vector<Vector3S> getForce(int eIdx);
  std::vector<Vector3S> getForce();

  MatrixXS getStiffness(int eIdx);

  ///@param trig if true, return only the upper triangle of the symmetric matrix.
  void getStiffnessSparse(std::vector<cfgScalar> &val, bool trig = false, bool constrained=false, bool iFixedRigid=false);
  void getStiffnessSparse(std::vector<cfgScalar> &val, bool trig, bool constrained, bool iFixedTranslation, bool iFixedRotation, bool iPeriodic);

  Eigen::SparseMatrix<cfgScalar> getStiffnessSparse(bool trig = false, bool constrained=false, bool iFixedRigid=false);

  ///@brief precompute lumped maasses.
  void computeMass();

  ///@param I row offsets. I.size() = matrix size + 1. I[size()-1]=number of non-zeros.
  ///@param J column indices.
  void stiffnessPattern(std::vector<int> & I, std::vector<int> & J, bool trig = false, bool iFixedRigid=false);
  void stiffnessPattern(std::vector<int> & I, std::vector<int> & J, bool trig, bool iFixedTranslation, bool iFixedRotation, bool iPeriodic);

  MatrixXS getStiffness();
  
  cfgScalar eleSize();

  virtual ~ElementMesh();


  //for rendering 
  cfgScalar forceDrawingScale;

//private:
  void fixTranslation(Eigen::SparseMatrix<cfgScalar> & K, bool iTriangular, ElementMesh * mesh);
  void fixRotation(Eigen::SparseMatrix<cfgScalar> & K, bool iTriangular, ElementMesh * mesh);
  void enforcePeriodicity(Eigen::SparseMatrix<cfgScalar> & K, bool iTriangular, ElementMesh * mesh);
};


void getEleX(int ii, const ElementMesh * m, std::vector<Vector3S> &x);

void setEleX(int ii, ElementMesh * m, const std::vector<Vector3S> &x);

#endif
