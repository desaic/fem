#include "Element2D.h"
#include "femError.hpp"

Element2D::Element2D(int _n):n(_n)
{}

Element2D::Element2D(const Element2D & e) : n(e.getNodeIndices())
{}
  
Element2D::Element2D(const std::vector<int> &inodes)
{
  n = inodes;
}

cfgScalar Element2D::getVol(const std::vector<Vector2S> & X)
{
  return 0;
}

Vector2S Element2D::getDisp(const Vector2S & p, const std::vector<Vector2S> & X,
    const std::vector<Vector2S>x)
{
  std::vector<cfgScalar> w = shapeFun(p);
  Vector2S u(0,0);
  for(unsigned int ii = 0; ii<w.size(); ii++){
    u += w[ii]*(x[ii] - X[ii]);
  }
  return u;
}

Matrix2S
Element2D::defGrad(Vector2S p, const std::vector<Vector2S> & X,
  const std::vector<Vector2S> & x) const
{
  Matrix2S F=Matrix2S::Identity();
  for(int ii = 0; ii<nV(); ii++){
    int vi = at(ii);
    Vector2S gradN = shapeFunGrad(ii,p,X);
    //outer product
    F += (x[vi] - X[vi])*gradN.transpose();
  }
  //fem_error = FEM_OK;
  cfgScalar det = F(0,0)*F(1,1)-F(0,1)*F(1,0);
  //if(F.determinant()<0)
  if (det<0)
  {
    fem_error = FEM_ERROR_INVERT;
  }
  return F;
}

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
