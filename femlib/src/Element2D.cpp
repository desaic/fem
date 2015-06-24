#include "Element2D.h"
#include "femError.hpp"

float Element2D::getVol(const std::vector<Vector2f> & X)
{
  return 0;
}

Vector2f Element2D::getDisp(const Vector2f & p, const std::vector<Vector2f> & X,
    const std::vector<Vector2f>x)
{
  std::vector<float> w = shapeFun(p);
  Vector2f u(0,0);
  for(unsigned int ii = 0; ii<w.size(); ii++){
    u += w[ii]*(x[ii] - X[ii]);
  }
  return u;
}

Matrix2f
Element2D::defGrad(Vector2f p, const std::vector<Vector2f> & X,
  const std::vector<Vector2f> & x) const
{
  Matrix2f F=Matrix2f::identity();
  for(int ii = 0; ii<nV(); ii++){
    int vi = at(ii);
    Vector2f gradN = shapeFunGrad(ii,p,X);
    //outer product
    F += outerProd((x[vi] - X[vi]) , gradN);
  }
  //fem_error = FEM_OK;
  if(F.determinant()<0){
    fem_error = FEM_ERROR_INVERT;
  }
  return F;
}

Element2D::Element2D(int _n):n(_n)
{}

Element2D::Element2D(const Element2D & e) : n(e.getNodeIndices())
{}
  

std::vector<std::array<int,2> >
Element2D::getEdges()
{
  std::vector<std::array<int,2> >edges;
  return edges;
}

int Element2D::nV()const
{
  return (int)n.size();
}

int Element2D::nF()const
{
	return 0;
}

int &
Element2D::operator[](int idx)
{
  return n[idx];
}

int
Element2D::operator[](int idx)const
{
  return n[idx];
}

int Element2D::at(int idx)const
{
  return n[idx];
}

void Element2D::resize(int size)
{
  n.resize(size);
}
