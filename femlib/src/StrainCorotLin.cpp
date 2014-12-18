#include "StrainCorotLin.hpp"
#include <Eigen\dense>

Eigen::Matrix3f toEigen(const Matrix3f & F)
{
  Eigen::Matrix3f m;
  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj < 3; jj++){
      m(ii, jj) = F(ii, jj);
    }
  }
  return m;
}

Matrix3f toVecmath(const Eigen::Matrix3f & F)
{
  Matrix3f m;
  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj < 3; jj++){
      m(ii, jj) = F(ii, jj);
    }
  }
  return m;
}

StrainCorotLin::StrainCorotLin()
{
  param.resize(2);
  param[0] = 1;
  param[1] = 10;
}

float StrainCorotLin::getEnergy(const Matrix3f & F)
{
  Eigen::Matrix3f m = toEigen(F);
  Eigen::JacobiSVD<Eigen::Matrix3f> svd(m);
  Eigen::Vector3f Sigma = svd.singularValues() - Eigen::Vector3f(1,1,1);
  float mu = param[0];
  float lambda = param[1];
  Eigen::Matrix3f I = Eigen::Matrix3f::Identity();
  float t = Sigma[0] + Sigma[1] + Sigma[2];
  float Psi = mu*(Sigma[0] * Sigma[0] + Sigma[1] * Sigma[1] + Sigma[2] * Sigma[2]) + 0.5f * lambda * t * t;
  return Psi;
}

Matrix3f StrainCorotLin::getPK1(const Matrix3f & F)
{
  Eigen::Matrix3f m = toEigen(F);
  Eigen::JacobiSVD<Eigen::Matrix3f> svd(m, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Eigen::Matrix3f U = svd.matrixU();
  Eigen::Matrix3f V = svd.matrixV();
  Eigen::Matrix3f R = U * V.transpose();
  float mu = param[0];
  float lambda = param[1];
  Eigen::Matrix3f I = Eigen::Matrix3f::Identity();
  Eigen::Matrix3f P = 2*mu*(m-R) + lambda * (R.transpose() * m - I).trace() * R;

  Matrix3f ret = toVecmath(P);
  return ret;
}

Matrix3f StrainCorotLin::getdPdx(const Matrix3f & F, const Matrix3f & dF)
{
  return Matrix3f();
}