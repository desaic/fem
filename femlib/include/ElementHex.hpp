#ifndef ELEMENTHEX_HPP
#define ELEMENTHEX_HPP
#include "Element.hpp"
class ElementHex:public Element
{
public:
  ElementHex();

  std::vector<float> ShapeFun(const Vector3f & p)const;

  Vector3f ShapeFunGrad(int ii, const Vector3f & xx,
    const std::vector<Vector3f> & X) const;

  std::vector<std::array<int,2> > GetEdges();
};
#endif