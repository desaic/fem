#include "StrainEneNeo.hpp"
using namespace Eigen;
StrainEneNeo::StrainEneNeo()
{
  param.resize(2);
  param[0] =1;
  param[1] =10;
}

Matrix3f StrainEneNeo::getStrainTensor(const Matrix3f & F)
{
  return F.transpose() * F;
}

Eigen::MatrixXf StrainEneNeo::EMatrix()
{
  Eigen::MatrixXf E = Eigen::MatrixXf::Zero(6, 6);
  float G = param[0];
  float l = param[1];
  float e = G*(3 * l + 2 * G) / (l + G);
  float nu = l / (2 * (l + G));
  float c = e / ((1 + nu)*(1 - 2 * nu));
  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj < 3; jj++){
      E(ii, jj) = nu;
    }
    E(ii, ii) = 1 - nu;
    E(3 + ii, 3 + ii) = 0.5f - nu;
  }
  return c*E;
}

float StrainEneNeo::getEnergy(const Matrix3f & F)
{
  float I1 = (F.transpose()*F).trace();
  float JJ = std::log(F.determinant());
  float mu = param[0],lambda=param[1];
  float Psi = (mu/2) * (I1-3) - mu*JJ + (lambda/2)*JJ*JJ;
  return (float)Psi;
}

Matrix3f StrainEneNeo::getPK1(const Matrix3f & F)
{
  float JJ = std::log(F.determinant());
  Matrix3f Finv = F.inverse();
  Finv.transpose();
  float mu = param[0],lambda=param[1];
  Matrix3f PP = mu*(F-Finv) + lambda*JJ*Finv;
  return PP;
}

Matrix3f StrainEneNeo::getdPdx(const Matrix3f & F,const Matrix3f & dF)
{
  Matrix3f dP = Matrix3f();
  float JJ = std::log(F.determinant());
  Matrix3f Finv = F.inverse();
  Matrix3f FinvT = Finv.transpose();
  float mu = param[0],lambda=param[1];
  dP = mu*dF;
  float c1 = mu-lambda * JJ;
  dP += c1 * FinvT*dF.transpose()*FinvT;
  dP += lambda*(Finv*dF).trace()*FinvT;
  return dP;
}
