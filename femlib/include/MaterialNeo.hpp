#ifndef MATERIALNEO_HPP
#define MATERIALNEO_HPP

#include "MaterialQuad.hpp"

#include <vecmath.h>

class MaterialNeo:public MaterialQuad
{
public:
  MaterialNeo();
  float getEnergy(const Matrix3f & F);
  Matrix3f getPK1(const Matrix3f & F);
  Matrix3f getdPdx(const Matrix3f & F,const Matrix3f & dF);
  
};
#endif