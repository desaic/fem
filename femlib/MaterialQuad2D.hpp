#ifndef MATERIALQUAD2D_HPP
#define MATERIALQUAD2D_HPP
#include "vecmath.h"
#include "Material2D.hpp"

class StrainEne;
class Quadrature;

///@brief abstract material model that uses quadrature
class MaterialQuad2D:public Material2D
{
public:
  MaterialQuad2D(StrainEne * ene = 0, Quadrature * _q = 0);

  MaterialQuad2D(const std::vector<StrainEne *> &ene, Quadrature * _q = 0);
  
  ///@brief precompute gradN.
  void init(ElementMesh2D * mesh);

  float getEnergy(Element2D* ele, ElementMesh2D * mesh);
  std::vector<Vector2f> getForce(Element2D* ele, ElementMesh2D * mesh);
  
  MatrixXf getStiffness(Element2D* ele, ElementMesh2D * mesh);

  std::vector<StrainEne*> e;
  const Quadrature * q;

  //pre-computed gradN for a unit-size element.
  //gradN[ii][q] is shape function gradient with respect to vertex ii and quadrature point q.
  //does not work for meshes with different element shape functions.
  std::vector<std::vector<Vector2> > gradN;

private:
  Eigen::MatrixXf stiffness(int qi, const MaterialQuad2D * mat, Element2D* ele, ElementMesh2D * mesh)
};

#endif