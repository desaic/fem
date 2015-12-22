#ifndef MATERIALQUAD_HPP
#define MATERIALQUAD_HPP
#include "Material.hpp"
#include <Eigen/Dense>
class StrainEne;
class Quadrature;

///@brief abstract material model that uses quadrature
class MaterialQuad:public Material
{
public:
  MaterialQuad(StrainEne * ene = 0, Quadrature * _q = 0);

  MaterialQuad(const std::vector<StrainEne *> &ene, Quadrature * _q = 0);
  
  ///@brief precompute gradN.
  void init(ElementMesh * mesh);

  float getEnergy(Element* ele, ElementMesh * mesh);
  std::vector<Eigen::Vector3f> getForce(Element* ele, ElementMesh * mesh);
  
  Eigen::MatrixXf getStiffness(Element* ele, ElementMesh * mesh);

  virtual std::vector<Eigen::MatrixXf> getElasticityTensors();
  virtual std::vector<Eigen::Matrix3f> getStrainTensors(Element* ele, ElementMesh * mesh, const std::vector<Eigen::Vector3f> &ix);

  std::vector<StrainEne*> e;
  const Quadrature * q;

  //pre-computed gradN for a unit-size element.
  //gradN[ii][q] is shape function gradient with respect to vertex ii and quadrature point q.
  //does not work for meshes with different element shape functions.
  std::vector<std::vector<Eigen::Vector3f> > gradN;
};

#endif