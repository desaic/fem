#ifndef MATERIAL_HPP
#define MATERIAL_HPP
#include <vector>
#include <Eigen/Dense>
#include "cfgDefs.h"

class Element;
class ElementMesh;

class Material
{
public:
  Material();
  virtual cfgScalar getEnergy(Element* ele, ElementMesh * mesh)=0;
  virtual std::vector<Vector3S> getForce(Element* ele, ElementMesh * mesh) = 0;
  virtual MatrixXS getStiffness(Element* ele, ElementMesh * mesh);
  virtual ~Material();

  virtual void init(ElementMesh * mesh){};
  virtual std::vector<MatrixXS> getElasticityTensors()=0;
  virtual std::vector<Matrix3S> getStrainTensors(Element* ele, ElementMesh* mesh, const std::vector<Vector3S> &ix) = 0;
  
  std::vector<cfgScalar> param;
};

#endif
