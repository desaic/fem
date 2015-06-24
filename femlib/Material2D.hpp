#ifndef MATERIAL2D_HPP
#define MATERIAL2D_HPP
#include "vecmath.h"
#include <vector>
#include "MatrixX.hpp"

class Element2D;
class ElementMesh2D;

class Material2D
{
public:
  Material2D();
  virtual float getEnergy(Element2D* ele, ElementMesh2D * mesh)=0;
  virtual std::vector<Vector2f> getForce(Element2D* ele, ElementMesh2D * mesh)=0;
  virtual MatrixXf getStiffness(Element2D* ele, ElementMesh2D * mesh);
  virtual ~Material2D();

  virtual void init(ElementMesh2D * mesh){};
  
  std::vector<float> param;
};

#endif