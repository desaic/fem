#include "StrainEneNeo.hpp"
using namespace Eigen;
StrainEneNeo::StrainEneNeo()
{
  param.resize(2);
  param[0] =1;
  param[1] =10;
}

Matrix3S StrainEneNeo::getStrainTensor(const Matrix3S & F)
{
  return F.transpose() * F;
}

MatrixXS StrainEneNeo::EMatrix()
{
  MatrixXS E = MatrixXS::Zero(6, 6);
  cfgScalar G = param[0];
  cfgScalar l = param[1];
  cfgScalar e = G*(3 * l + 2 * G) / (l + G);
  cfgScalar nu = l / (2 * (l + G));
  cfgScalar c = e / ((1 + nu)*(1 - 2 * nu));
  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj < 3; jj++){
      E(ii, jj) = nu;
    }
    E(ii, ii) = 1 - nu;
    E(3 + ii, 3 + ii) = 0.5f - nu;
  }
  return c*E;
}

cfgScalar StrainEneNeo::getEnergy(const Matrix3S & F)
{
  cfgScalar I1 = (F.transpose()*F).trace();
  cfgScalar JJ = std::log(F.determinant());
  cfgScalar mu = param[0],lambda=param[1];
  cfgScalar Psi = (mu/2) * (I1-3) - mu*JJ + (lambda/2)*JJ*JJ;
  return (cfgScalar)Psi;
}

Matrix3S StrainEneNeo::getPK1(const Matrix3S & F)
{
  cfgScalar JJ = std::log(F.determinant());
  Matrix3S Finv = F.inverse().transpose();
  cfgScalar mu = param[0],lambda=param[1];
  Matrix3S PP = mu*(F-Finv) + lambda*JJ*Finv;
  return PP;
}

Matrix3S StrainEneNeo::getdPdx(const Matrix3S & F,const Matrix3S & dF)
{
  Matrix3S dP = Matrix3S();
  cfgScalar JJ = std::log(F.determinant());
  Matrix3S Finv = F.inverse();
  Matrix3S FinvT = Finv.transpose();
  cfgScalar mu = param[0],lambda=param[1];
  dP = mu*dF;
  cfgScalar c1 = mu-lambda * JJ;
  dP += c1 * FinvT*dF.transpose()*FinvT;
  dP += lambda*(Finv*dF).trace()*FinvT;
  return dP;
}
