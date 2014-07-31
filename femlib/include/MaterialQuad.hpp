#ifndef MATERIALQUAD_HPP
#define MATERIALQUAD_HPP
#include "vecmath.h"
#include "Material.hpp"

class Quadrature;

///@brief abstract material model that uses quadrature
class MaterialQuad:public Material
{
public:
  MaterialQuad();
  float getEnergy(Element* ele, ElementMesh * mesh);
  std::vector<Vector3f> getForce(Element* ele, ElementMesh * mesh);
  MatrixXd getStiffness(Element* ele, ElementMesh * mesh);
  
  virtual float getEnergy(const Matrix3f & F)=0;
  virtual Matrix3f getPK1(const Matrix3f & F)=0;
  virtual Matrix3f getdPdx(const Matrix3f & F,const Matrix3f & dF)=0;
  virtual ~MaterialQuad();

  const Quadrature * q;
};

#endif