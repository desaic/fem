#ifndef MATERIAL_HPP
#define MATERIAL_HPP
#include "vecmath.h"
#include <vector>
class Element;
class ElementMesh;
struct MatrixXd;
class Material
{
public:
  Material();
  virtual float getEnergy(Element* ele, ElementMesh * mesh)=0;
  virtual std::vector<Vector3f> getForce(Element* ele, ElementMesh * mesh)=0;
  virtual MatrixXd getStiffness(Element* ele, ElementMesh * mesh);
  virtual ~Material();

  virtual void init(ElementMesh * mesh){};
  
  std::vector<float> param;
};

#endif