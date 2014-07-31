#include "Element.hpp"
float Element::getVol(const std::vector<Vector3f> & X)
{
  return 0;
}

Matrix3f
Element::defGrad(Vector3f p, const std::vector<Vector3f> & X,
  const std::vector<Vector3f> & x) const
{
  Matrix3f F=Matrix3f::identity();
  for(int ii = 0; ii<nV(); ii++){
    int vi = at(ii);
    Vector3f gradN = shapeFunGrad(ii,p,X);
    //outer product
    F += outerProd((x[vi] - X[vi]) , gradN);
  }
  return F;
}

Element::Element(int _n):n(_n)
{}

std::vector<std::array<int,2> >
Element::getEdges()
{
  std::vector<std::array<int,2> >edges;
  return edges;
}

int Element::nV()const
{
  return (int)n.size();
}

int &
Element::operator[](int idx)
{
  return n[idx];
}

int
Element::operator[](int idx)const
{
  return n[idx];
}

int Element::at(int idx)const
{
  return n[idx];
}

void Element::resize(int size)
{
  n.resize(size);
}
