#include "Material.hpp"
#include "MatrixXd.hpp"
Material::Material(){}
Material::~Material(){}
MatrixXd Material::getStiffness(Element* ele, ElementMesh * mesh)
{
  return MatrixXd();
}