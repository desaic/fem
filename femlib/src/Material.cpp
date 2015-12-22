#include "Material.hpp"
Material::Material(){}
Material::~Material(){}
Eigen::MatrixXf Material::getStiffness(Element* ele, ElementMesh * mesh)
{
  return Eigen::MatrixXf();
}