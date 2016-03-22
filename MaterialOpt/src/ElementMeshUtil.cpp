#include "ElementMeshUtil.hpp"
#include "Element.hpp"
#include "ElementMesh.hpp"
#include "cfgDefs.h"
#include "ArrayUtil.hpp"

static const int sw[8][3] =
{ { -1, -1, -1 },
{ -1, -1, 1 },
{ -1, 1, -1 },
{ -1, 1, 1 },
{ 1, -1, -1 },
{ 1, -1, 1 },
{ 1, 1, -1 },
{ 1, 1, 1 }
};

std::vector<int> topVerts(ElementMesh * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int ny = s[1];
  int nz = s[2];
  int topV[4] = { 2, 3, 6, 7 };
  for (int ii = 0; ii<nx; ii++){
    for (int kk = 0; kk < nz; kk++){
      int ei = gridToLinearIdx(ii, ny - 1, kk, s);
      for (int jj = 0; jj < 4; jj++){
        int vi = em->e[ei]->at(topV[jj]);
        v.push_back(vi);
      }
    }
  }
  return v;
}

std::vector<int> botVerts(ElementMesh * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int nz = s[2];
  int botV[4] = { 0, 1, 4, 5};
  for (int ii = 0; ii<nx; ii++){
    for (int kk = 0; kk < nz; kk++){
      int ei = gridToLinearIdx(ii,0,kk,s);
      for (int jj = 0; jj<4; jj++){
        int vi = em->e[ei]->at(botV[jj]);
        v.push_back(vi);
      }
    }
  }
  return v;
}

std::vector<int> leftVerts(ElementMesh * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int ny = s[1];
  int nz = s[2];
  int leftV[4] = { 0, 1 ,2 ,3};
  for (int ii = 0; ii<ny; ii++){
    for (int kk = 0; kk < nz; kk++){
      int ei = gridToLinearIdx(0,ii,kk, s);
      for (int jj = 0; jj < 4; jj++){
        int vi = em->e[ei]->at(leftV[jj]);
        v.push_back(vi);
      }
    }
  }
  return v;
}

std::vector<int> rightVerts(ElementMesh * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int ny = s[1];
  int nz = s[2];
  int rightV[4] = { 4, 5, 6, 7 };
  for (int ii = 0; ii<ny; ii++){
    for (int kk = 0; kk < nz; kk++){
      int ei = gridToLinearIdx(nx-1,ii,kk,s);
      for (int jj = 0; jj < 4; jj++){
        int vi = em->e[ei]->at(rightV[jj]);
        v.push_back(vi);
      }
    }
  }
  return v;
}

std::vector<int> frontVerts(ElementMesh * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int ny = s[1];
  int nz = s[2];
  int rightV[4] = { 1, 3, 5, 7 };
  for (int ii = 0; ii<nx; ii++){
    for (int jj = 0; jj < ny; jj++){
      int ei = gridToLinearIdx(ii, jj, nz-1, s);
      for (int jj = 0; jj < 4; jj++){
        int vi = em->e[ei]->at(rightV[jj]);
        v.push_back(vi);
      }
    }
  }
  return v;
}

std::vector<int> backVerts(ElementMesh * em, const std::vector<int> & s)
{
  std::vector<int> v;
  int nx = s[0];
  int ny = s[1];
  int rightV[4] = { 0, 2, 4, 6 };
  for (int ii = 0; ii<nx; ii++){
    for (int jj = 0; jj < ny; jj++){
      int ei = gridToLinearIdx(ii, jj, 0, s);
      for (int jj = 0; jj < 4; jj++){
        int vi = em->e[ei]->at(rightV[jj]);
        v.push_back(vi);
      }
    }
  }
  return v;
}


void stretchX(ElementMesh * em, const Eigen::Vector3d & ff, const std::vector<int> & s, std::vector<double> & externalForce)
{
  std::vector<int> leftv, rightv;
  rightv = rightVerts(em, s);
  Eigen::Vector3d fv = ff / (double)rightv.size();
  addVector3d(externalForce, fv, rightv);
  leftv = leftVerts(em, s);
  addVector3d(externalForce, -fv, leftv);
}

void stretchY(ElementMesh * em, const Eigen::Vector3d & ff, const std::vector<int> & s, std::vector<double> & externalForce)
{
  std::vector<int> topv, botv;
  topv = topVerts(em, s);
  Eigen::Vector3d fv = ff / (double)topv.size();
  addVector3d(externalForce, fv, topv);
  botv = botVerts(em, s);
  addVector3d(externalForce, -fv, botv);
}

void stretchZ(ElementMesh * em, const Eigen::Vector3d & ff, const std::vector<int>& s, std::vector<double> & externalForce)
{
  std::vector<int> frontv, backv;
  frontv = frontVerts(em, s);
  Eigen::Vector3d fv = ff / (double)frontv.size();
  addVector3d(externalForce, fv, frontv);
  backv = backVerts(em, s);
  addVector3d(externalForce, -fv, backv);
}

void shearXY(ElementMesh * em, double ff,
    const std::vector<int>& s, std::vector<double> & fe)
{
  Eigen::Vector3d force = ff * Eigen::Vector3d(1, 0, 0);
  //apply horizontal force on top and bottom faces
  stretchY(em, force, s, fe);
  force = ff * Eigen::Vector3d(0, 1, 0);
  stretchX(em, force, s, fe);
}

void shearYZ(ElementMesh * em, double ff,
  const std::vector<int>& s, std::vector<double> & fe)
{
  Eigen::Vector3d force = ff * Eigen::Vector3d(0, 1, 0);
  //apply horizontal force on top and bottom faces
  stretchZ(em, force, s, fe);
  force = ff * Eigen::Vector3d(0, 0, 1);
  stretchY(em, force, s, fe);
}

void shearXZ(ElementMesh * em, double ff, 
  const std::vector<int>& s, std::vector<double> & fe)
{
  Eigen::Vector3d force = ff * Eigen::Vector3d(1, 0, 0);
  //apply horizontal force on top and bottom faces
  stretchZ(em, force, s, fe);
  force = ff * Eigen::Vector3d(0, 0, 1);
  stretchX(em, force, s, fe);
}

std::vector<int> cornerVerts(ElementMesh * em, const std::vector<int> & gridSize)
{
  std::vector<int> corner(8, 0);
  int x = gridSize[0] - 1;
  int y = gridSize[1] - 1;
  int z = gridSize[2] - 1;
  corner[0] = em->e[gridToLinearIdx(0, 0, 0, gridSize)]->at(0);
  corner[1] = em->e[gridToLinearIdx(0, 0, z, gridSize)]->at(1);
  corner[2] = em->e[gridToLinearIdx(0, y, 0, gridSize)]->at(2);
  corner[3] = em->e[gridToLinearIdx(0, y, z, gridSize)]->at(3);
  corner[4] = em->e[gridToLinearIdx(x, 0, 0, gridSize)]->at(4);
  corner[5] = em->e[gridToLinearIdx(x, 0, z, gridSize)]->at(5);
  corner[6] = em->e[gridToLinearIdx(x, y, 0, gridSize)]->at(6);
  corner[7] = em->e[gridToLinearIdx(x, y, z, gridSize)]->at(7);
  return corner;
}

//assuming element size 1.
Vector3d shapeFunGrad(int ii, const Vector3d & xx)
{
  Vector3d grad;
  grad[0] = sw[ii][0] * (1 + sw[ii][1] * xx[1]) * (1 + sw[ii][2] * xx[2]);
  grad[1] = sw[ii][1] * (1 + sw[ii][0] * xx[0]) * (1 + sw[ii][2] * xx[2]);
  grad[2] = sw[ii][2] * (1 + sw[ii][0] * xx[0]) * (1 + sw[ii][1] * xx[1]);
  return 0.25*grad;
}

Eigen::MatrixXd BMatrix(const Vector3d & xx, const Eigen::Vector3d & size)
{
  int nV = 8;
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(6, 3 * nV);
  for (int ii = 0; ii < nV; ii++){
    int col = 3 * ii;
    Vector3d dN = shapeFunGrad(ii, xx).cwiseQuotient(size);
    B(0, col) = dN[0];
    B(1, col + 1) = dN[1];
    B(2, col + 2) = dN[2];

    B(3, col) = dN[1];
    B(3, col + 1) = dN[0];

    B(4, col + 1) = dN[2];
    B(4, col + 2) = dN[1];

    B(5, col) = dN[2];
    B(5, col + 2) = dN[0];
  }
  return B;
}

//returns 6 element vector measured at xi.
Eigen::VectorXd hexStrain(const Eigen::VectorXd & u, const Eigen::VectorXd & X,
  const Eigen::Vector3d & xi)
{
  int dim = 3;
  //size of cube element
  int nV = X.rows()/3;
  Eigen::Vector3d esize = X.segment<3>(3*(nV-1)) - X.segment<3>(0);
  Eigen::MatrixXd B=BMatrix(xi, esize);
  Eigen::VectorXd strain = B*u;
  return strain;
}
