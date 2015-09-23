#ifndef MATERIAL2D_HPP
#define MATERIAL2D_HPP
#include "vecmath.h"
#include <vector>
#include "MatrixX.hpp"
#include <Eigen\Dense>
#include "cfgDefs.h"

class Element2D;
class ElementMesh2D;

class Material2D
{
public:
  Material2D();
  virtual ~Material2D();

  virtual void init(ElementMesh2D * mesh){};

  virtual cfgScalar getEnergy(Element2D* ele, ElementMesh2D * mesh)=0;
  virtual std::vector<Vector2S> getForce(Element2D* ele, ElementMesh2D * mesh)=0;
  virtual MatrixXS getStiffness(Element2D* ele, ElementMesh2D * mesh);
  virtual std::vector<Matrix2S> getStressTensors(Element2D* ele, ElementMesh2D * mesh)=0;
  virtual std::vector<Matrix2S> getStrainTensors(Element2D* ele, ElementMesh2D * mesh) = 0;
  virtual std::vector<Matrix2S> getStrainTensors(Element2D* ele, ElementMesh2D * mesh, const std::vector<Vector2S> &ix) = 0;

  virtual std::vector<MatrixXS> getElasticityTensors()=0;

  std::vector<cfgScalar> param;
};

#endif