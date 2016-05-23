#ifndef ELEMENTHEX_HPP
#define ELEMENTHEX_HPP
#include "Element.hpp"

#include <Eigen/Dense>
extern const int facen[][3];
extern const int faces[][4];
class ElementHex:public Element
{
public:
  ElementHex();
  ElementHex(const ElementHex & e);
  ElementHex(const std::vector<int> &iNodes);

  int nV() const override{ return 8; }

  int nF() const override{ return 6; }


  ///@brief natural coordinate for a point in reference space
  Vector3S natCoord(const Vector3S & p, const std::vector<Vector3S> & X);
  
  std::vector<cfgScalar> shapeFun(const Vector3S & p)const;

  Vector3S shapeFunGrad(int ii, const Vector3S & xx,
    const std::vector<Vector3S> & X) const;

  std::vector<std::array<int,2> > getEdges();
  cfgScalar getVol(const std::vector<Vector3S> & X);

  Vector3S facePt(int fi, const Vector3S & x);
  ///@param fi face index.
  MatrixXS NMatrix(int fi);
  ///@brief shape function matrix.
  MatrixXS HMatrix(const Vector3S & xx);
  //specific for linear strain (F+F^T).
  MatrixXS BMatrix(const Vector3S & xx, const std::vector<Vector3S>X);
};
#endif
