#include "MaterialNeo.hpp"
#include "MatrixXd.hpp"
#include "Quadrature.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"

MaterialNeo::MaterialNeo()
{
  param.resize(2);
  param[0] =1;
  param[1] =10;
  q=&(Quadrature::Gauss2);
}

float MaterialNeo::getEnergy(const Matrix3f & F)
{
  float I1 = (F.transposed()*F).trace();
  float JJ = std::log(F.determinant());
  float mu = param[0],lambda=param[1];
  float Psi = (mu/2) * (I1-3) - mu*JJ + (lambda/2)*JJ*JJ;
  return (float)Psi;
}

float MaterialNeo::getEnergy(Element* ele, ElementMesh * mesh)
{
  float energy = 0;
  for(int ii = 0; ii<q->x.size();ii++){
    Matrix3f F = ele->defGrad(q->x[ii], mesh->X, mesh->x);
    energy += q->w[ii] * getEnergy(F);
  }
  return ele->getVol(mesh->X) * energy;
}

std::vector<Vector3f> MaterialNeo::getForce(Element* ele, ElementMesh * mesh)
{

}

Matrix3f MaterialNeo::getPK1(const Matrix3f & F)
{
  float JJ = std::log(F.determinant());
  Matrix3f Finv = F.inverse();
  Finv.transpose();
  float mu = param[0],lambda=param[1];
  Matrix3f PP = mu*(F-Finv) + lambda*JJ*Finv;
  return PP;
}

Matrix3f MaterialNeo::getdPdx(const Matrix3f & F,const Matrix3f & dF)
{
}

MatrixXd MaterialNeo::getStiffness(Element* ele, ElementMesh * mesh)
{

}
