#include "UnitTests.hpp"
#include "ElementRegGrid.hpp"
#include "ElementHier.hpp"
#include <iostream>
#include <vector>

///@param eidx fine element index
///@param p natural coordinate in fine element
///@param m all elements are assumed to be ElementHier
Vector3f getDisp(int eidx, std::vector<ElementMesh *> & meshes,  const Vector3f & p)
{
  //coordinate of point in reference space
  Vector3f Xp(0, 0, 0);
  ElementHier * fine = (ElementHier*)meshes[0]->e[eidx];
  std::vector<float> N = fine->shapeFun(p);
  for (int ii = 0; ii < 8; ii++){
    Xp += meshes[0]->X[fine->at(ii)] * N[ii];
  }

  Vector3f u(0, 0, 0);
  int ei = eidx;
  for (unsigned int ii = 0; ii < meshes.size() - 1; ii++){
    ElementMesh * m = meshes[ii];
    ElementHier * e = (ElementHier*)m->e[ei];
    //natural coordiante of p in current level of element;
    Vector3f Xn = e->natCoord(Xp, m->X);
    N = e->shapeFun(Xn);
    //displacement of one level
    Vector3f ul(0, 0, 0);
    for (int jj = 0; jj < 8; jj++){
      ul += N[jj]*(m->X[e->at(jj)] - m->x[e->at(jj)]);
    }
    Matrix3f F = e->defGrad(Xp, m->X, m->x);
    u = ul + F*u;
    ei = e->parent;
  }
  return u;
}

Matrix3f defGrad(int eidx, std::vector<ElementMesh *> & meshes, const Vector3f & p)
{
  //coordinate of point in reference space
  Vector3f Xp(0, 0, 0);
  ElementHier * fine = (ElementHier*)meshes[0]->e[eidx];
  std::vector<float> N = fine->shapeFun(p);
  for (int ii = 0; ii < 8; ii++){
    Xp += meshes[0]->X[fine->at(ii)] * N[ii];
  }

  Matrix3f F = Matrix3f::identity();
  int ei = eidx;
  for (unsigned int ii = 0; ii < meshes.size() - 1; ii++){
    ElementMesh * m = meshes[ii];
    ElementHier * e = (ElementHier*)m->e[ei];
    //natural coordiante of p in current level of element;
    Vector3f Xn = e->natCoord(Xp, m->X);
    Matrix3f Fl = e->defGrad(Xn, m->X, m->x);
    F = Fl*F;
  }
  return F;
}

void testCoarseDefGrad()
{
  ElementRegGrid *grid = new ElementRegGrid(2, 4, 2);
  replaceElementHex(grid);
  ElementMesh * cm = coarsen(grid);
  for (unsigned int ii = 0; ii < grid->e.size(); ii++){
    ElementHier * fine = (ElementHier*)grid->e[ii];
    std::cout << fine->parent<<"\n";
    for (unsigned int jj = 0; jj < fine->Xn.size(); jj++){
      std::cout << fine->Xn[jj][0] << " " << fine->Xn[jj][1] << " " << fine->Xn[jj][2] << "\n";
    }
  }

  for (unsigned int ii = 0; ii < cm->X.size(); ii++){
    std::cout << cm->X[ii][0] << " " << cm->X[ii][1] << " " << cm->X[ii][2] << "\n";
  }

  std::vector<ElementMesh * > m;
  m.push_back(grid);
  m.push_back(cm);
  Vector3f p(-0.57f, 0.57f, -0.57f);
  Vector3f u = getDisp(0, m, p);
  std::cout << u[0] << " " << u[1] << " " << u[2] << "\n";
}

void ElementCoarseTest()
{
  ElementRegGrid grid(2, 2, 2);
  ElementMesh* cm;
  cm = coarsen(&grid);

}
