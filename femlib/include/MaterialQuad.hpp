#ifndef MATERIALQUAD_HPP
#define MATERIALQUAD_HPP
#include "Material.hpp"
#include <Eigen/Dense>
#include "cfgDefs.h"
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

  cfgScalar getEnergy(Element* ele, ElementMesh * mesh);
  std::vector<Vector3S> getForce(Element* ele, ElementMesh * mesh);
  
  MatrixXS getStiffness(Element* ele, ElementMesh * mesh);

  virtual std::vector<MatrixXS> getElasticityTensors();
  virtual std::vector<Matrix3S> getStrainTensors(Element* ele, ElementMesh * mesh, const std::vector<Vector3S> &ix);

  std::vector<StrainEne*> e;
  const Quadrature * q;

  //pre-computed gradN for a unit-size element.
  //gradN[ii][q] is shape function gradient with respect to vertex ii and quadrature point q.
  //does not work for meshes with different element shape functions.
  std::vector<std::vector<Vector3S> > gradN;
};

#endif