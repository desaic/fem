#include "Material2D.h"
#include "MatrixX.hpp"
Material2D::Material2D(){}
Material2D::~Material2D(){}
MatrixXf Material2D::getStiffness(Element2D* ele, ElementMesh2D * mesh)
{
  return MatrixXf();
}


