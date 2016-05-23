#include "Element.hpp"
#include "femError.hpp"
using namespace Eigen;
cfgScalar Element::getVol(const std::vector<Vector3S> & X)
{
  return 0;
}

Vector3S Element::getDisp(const Vector3S & p, const std::vector<Vector3S> & X,
    const std::vector<Vector3S>x)
{
  std::vector<cfgScalar> w = shapeFun(p);
  Vector3S u(0,0,0);
  for(unsigned int ii = 0; ii<w.size(); ii++){
    u += w[ii]*(x[ii] - X[ii]);
  }
  return u;
}

Matrix3S
Element::defGrad(Vector3S p, const std::vector<Vector3S> & X,
  const std::vector<Vector3S> & x) const
{
  Matrix3S F=Matrix3S::Identity();
  for(int ii = 0; ii<nV(); ii++){
    int vi = at(ii);
    Vector3S gradN = shapeFunGrad(ii,p,X);
    //outer product
    F += (x[vi] - X[vi]) * (gradN.transpose());
  }
  //fem_error = FEM_OK;
  if(F.determinant()<0){
    fem_error = FEM_ERROR_INVERT;
  }
  return F;
}

Element::Element(int _n):n(_n)
{}

Element::Element(const Element & e) : n(e.getNodeIndices())
{}
  
Element::Element(const std::vector<int> &inodes)
{
  n = inodes;
}

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

int Element::nF()const
{
	return 0;
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
