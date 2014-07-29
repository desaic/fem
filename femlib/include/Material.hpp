#ifndef MATERIAL_HPP
#define MATERIAL_HPP
#include "vecmath.h"
#include <vector>
class Element;
class ElementMesh;
class Material
{
public:
  Material();
  virtual float getEnergy(Element* ele, ElementMesh * mesh)=0;
  virtual Matrix3f getForce(Element* ele, ElementMesh * mesh)=0;
  virtual float getPK1(Element* ele, ElementMesh * mesh)=0;
  virtual ~Material();
};

#endif