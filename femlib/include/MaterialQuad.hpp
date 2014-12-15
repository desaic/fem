#ifndef MATERIALQUAD_HPP
#define MATERIALQUAD_HPP
#include "vecmath.h"
#include "Material.hpp"

class StrainEne;
class Quadrature;

///@brief abstract material model that uses quadrature
class MaterialQuad:public Material
{
public:
  MaterialQuad(StrainEne * ene = 0, Quadrature * _q = 0);

  MaterialQuad(const std::vector<StrainEne *> &ene, Quadrature * _q = 0);
  
  float getEnergy(Element* ele, ElementMesh * mesh);
  std::vector<Vector3f> getForce(Element* ele, ElementMesh * mesh);
  
  MatrixXd getStiffness(Element* ele, ElementMesh * mesh);

  std::vector<StrainEne*> e;
  const Quadrature * q;
};

#endif