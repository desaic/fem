#ifndef MATERIALNEO_HPP
#define MATERIALNEO_HPP
#include "Material.hpp"
#include <vecmath.h>
class Quadrature;
class MaterialNeo:public Material
{
public:
  MaterialNeo();
  float getEnergy(Element* ele, ElementMesh * mesh);
  float getEnergy(const Matrix3f & F);
  std::vector<Vector3f> getForce(Element* ele, ElementMesh * mesh);
  Matrix3f getPK1(const Matrix3f & F);
  MatrixXd getStiffness(Element* ele, ElementMesh * mesh);
  Matrix3f getdPdx(const Matrix3f & F,const Matrix3f & dF);
  const Quadrature * q;
};
#endif