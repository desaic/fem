#ifndef MATERIAL_HPP
#define MATERIAL_HPP
#include <vector>
#include <Eigen/Dense>

class Element;
class ElementMesh;

class Material
{
public:
  Material();
  virtual float getEnergy(Element* ele, ElementMesh * mesh)=0;
  virtual std::vector<Eigen::Vector3f> getForce(Element* ele, ElementMesh * mesh) = 0;
  virtual Eigen::MatrixXf getStiffness(Element* ele, ElementMesh * mesh);
  virtual ~Material();

  virtual void init(ElementMesh * mesh){};
  virtual std::vector<Eigen::MatrixXf> getElasticityTensors()=0;
  virtual std::vector<Eigen::Matrix3f> getStrainTensors(Element* ele, ElementMesh* mesh, const std::vector<Eigen::Vector3f> &ix) = 0;
  
  std::vector<float> param;
};

#endif
