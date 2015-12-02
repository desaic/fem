#ifndef MATERIALQUAD2D_HPP
#define MATERIALQUAD2D_HPP
#include "vecmath.h"
#include "Material2D.h"
#include <Eigen/Dense>

class StrainEne2D;
class Quadrature2D;

///@brief abstract material model that uses quadrature
class MaterialQuad2D:public Material2D
{
public:
  MaterialQuad2D(StrainEne2D * ene = 0, Quadrature2D * _q = 0);

  MaterialQuad2D(const std::vector<StrainEne2D *> &ene, Quadrature2D * _q = 0);
  
  ///@brief precompute gradN.
  void init(ElementMesh2D * mesh);

  cfgScalar getEnergy(Element2D* ele, ElementMesh2D * mesh);
  std::vector<Vector2S> getForce(Element2D* ele, ElementMesh2D * mesh);
  
  MatrixXS getStiffness(Element2D* ele, ElementMesh2D * mesh);

  virtual std::vector<Matrix2S> getStressTensors(Element2D* ele, ElementMesh2D * mesh);
  virtual std::vector<Matrix2S> getStrainTensors(Element2D* ele, ElementMesh2D * mesh);

  virtual std::vector<Matrix2S> getStrainTensors(Element2D* ele, ElementMesh2D * mesh, const std::vector<Vector2S> &ix);

  virtual std::vector<MatrixXS> getElasticityTensors();

  std::vector<StrainEne2D*> e;
  const Quadrature2D * q;

  //pre-computed gradN for a unit-size element.
  //gradN[ii][q] is shape function gradient with respect to vertex ii and quadrature point q.
  //does not work for meshes with different element shape functions.
  std::vector<std::vector<Vector2S> > gradN;

private:
  MatrixXS stiffness(int qi, const MaterialQuad2D * mat, Element2D* ele, ElementMesh2D * mesh);
  bool m_planeStress;
};

#endif