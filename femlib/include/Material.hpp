#ifndef MATERIAL_HPP
#define MATERIAL_HPP
#include "vecmath.h"
#include <vector>
#include "MatrixX.hpp"
#include <Eigen\Dense>

class Element;
class ElementMesh;

class Material
{
public:
  Material();
  virtual float getEnergy(Element* ele, ElementMesh * mesh)=0;
  virtual std::vector<Vector3f> getForce(Element* ele, ElementMesh * mesh)=0;
  virtual MatrixXf getStiffness(Element* ele, ElementMesh * mesh);
  virtual ~Material();

  virtual void init(ElementMesh * mesh){};
  virtual std::vector<Eigen::MatrixXf> getElasticityTensors()=0;
  virtual std::vector<Matrix3f> getStrainTensors(Element* ele, ElementMesh* mesh, const std::vector<Vector3f> &ix)=0;
  
  std::vector<float> param;
};

#endif