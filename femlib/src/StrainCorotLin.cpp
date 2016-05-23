#include "StrainCorotLin.hpp"
#include <Eigen/Dense>
#include <vector>

#include <iostream>

std::vector<Matrix3S> initEp3();
//3D Levi-Civita symbol
std::vector<Matrix3S> Ep3=initEp3();
using namespace Eigen;
std::vector<Matrix3S> initEp3()
{
  std::vector<Matrix3S> Ep3(3, Matrix3S::Zero());
  Ep3[0](1, 2) = 1;
  Ep3[0](2, 1) = -1;

  Ep3[1](0, 2) = -1;
  Ep3[1](2, 0) = 1;

  Ep3[2](0, 1) = 1;
  Ep3[2](1, 0) = -1;

  return Ep3;
}

cfgScalar dot(Matrix3S m1, Matrix3S m2){
  cfgScalar prod = 0;
  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj < 3; jj++){
      prod += m1(ii, jj) * m2(ii, jj);
    }
  }
  return prod;
}

Matrix3S toEigen(const Matrix3S & F)
{
  Matrix3S m;
  for (int ii = 0; ii < 3; ii++){
    for (int jj = 0; jj < 3; jj++){
      m(ii, jj) = F(ii, jj);
    }
  }
  return m;
}

Matrix3S toVecmath(const Matrix3S & F)
{
  Matrix3S m;
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

cfgScalar StrainCorotLin::getEnergy(const Matrix3S & F)
{
  Matrix3S m = toEigen(F);
  Eigen::JacobiSVD<Matrix3S> svd(m);
  Vector3S Sigma = svd.singularValues() - Vector3S(1,1,1);
  cfgScalar mu = param[0];
  cfgScalar lambda = param[1];
  Matrix3S I = Matrix3S::Identity();
  cfgScalar t = Sigma[0] + Sigma[1] + Sigma[2];
  cfgScalar Psi = mu*(Sigma[0] * Sigma[0] + Sigma[1] * Sigma[1] + Sigma[2] * Sigma[2]) + 0.5f * lambda * t * t;
  return Psi;
}

Matrix3S StrainCorotLin::getPK1(const Matrix3S & F)
{
  Matrix3S m = toEigen(F);
  Eigen::JacobiSVD<Matrix3S> svd(m, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Matrix3S U = svd.matrixU();
  Matrix3S V = svd.matrixV();
  Matrix3S R = U * V.transpose();
  cfgScalar mu = param[0];
  cfgScalar lambda = param[1];
  Matrix3S I = Matrix3S::Identity();
  Matrix3S P = 2*mu*(m-R) + lambda * (R.transpose() * m - I).trace() * R;

  Matrix3S ret = toVecmath(P);
  return ret;
}

Matrix3S crossProdMat(const Vector3S & v)
{
  Matrix3S A = Matrix3S::Zero();
  A(0, 1) = -v[2];
  A(0, 2) = v[1];
  A(1, 0) = v[2];
  A(1, 2) = -v[0];
  A(2, 0) = -v[1];
  A(2, 1) = v[0];
  return A;
}

Matrix3S StrainCorotLin::getdPdx(const Matrix3S & m, const Matrix3S & dm)
{
  Matrix3S F  = toEigen(m);
  Matrix3S dF = toEigen(dm);
  Eigen::JacobiSVD<Matrix3S> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
  Matrix3S U = svd.matrixU();
  Matrix3S V = svd.matrixV();
  Matrix3S R = U * V.transpose();
  Matrix3S W = R.transpose()*dF;
  Matrix3S Sigma = svd.singularValues().asDiagonal();
  Matrix3S S = V*Sigma*V.transpose();
  Vector3S w;
  Matrix3S I = Matrix3S::Identity();
  w[0] = W(1, 2) - W(2, 1);
  w[1] = W(2, 0) - W(0, 2);
  w[2] = W(0, 1) - W(1, 0);
  Matrix3S SI = (S.trace()*I - S).inverse();
  Matrix3S dR = -R*crossProdMat(SI*w);
  //debug dR computation
  //cfgScalar h = 0.01;
  //Matrix3S F1 = F + h*dF;
  //svd.compute(F1, Eigen::ComputeFullU | Eigen::ComputeFullV);
  //U = svd.matrixU();
  //V = svd.matrixV();
  //Matrix3S R1 = U * V.transpose();
  //R1 -= R;
  //R1 = (1/h)*R1;
  //std::cout << dR << "\n" << R1 << "\n\n";

  Matrix3S dP;
  cfgScalar mu = param[0];
  cfgScalar lambda = param[1];
  dP = 2 * mu*dF + lambda*W.trace()*R +(lambda*(S - I).trace() - 2 * mu)*dR;
  return toVecmath(dP);
}

MatrixXS StrainCorotLin::EMatrix()
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
