#ifndef ELEMENTHEX_HPP
#define ELEMENTHEX_HPP
#include "Element.hpp"

#include <Eigen\Dense>

class ElementHex:public Element
{
public:
  ElementHex();

  std::vector<float> shapeFun(const Vector3f & p)const;

  Vector3f shapeFunGrad(int ii, const Vector3f & xx,
    const std::vector<Vector3f> & X) const;

  std::vector<std::array<int,2> > getEdges();
  float getVol(const std::vector<Vector3f> & X);

  //specific for linear strain (F+F^T)
  Eigen::MatrixXf BMatrix(const Vector3f & xx, const std::vector<Vector3f>X);
};
#endif