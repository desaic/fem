#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <array>
#include <vector>

#include "vecmath.h"

class Element
{
public:
  Element(int _n=0);
  ///@brief number of vertices
  int nV() const;
  int & operator[](int idx);
  int operator[](int idx)const;
  int at(int idx)const;
  void resize(int size);

  ///@param X Reference vertex positions.
  ///@param x global array of Deformed vertex positions.
  Matrix3f defGrad(Vector3f p,
    const std::vector<Vector3f> & x, const std::vector<Vector3f> & X) const;

  ///@param p natural coordinates.
  ///@return interpolation weight for each dof.
  virtual std::vector<float> shapeFun(const Vector3f & p) const=0;

  ///@brief ii index of basis function.
  ///@brief xx point in natural coordinates.
  ///@brief X global array of rest positions.
  virtual Vector3f shapeFunGrad(int ii, const Vector3f & xx,
    const std::vector<Vector3f> & X) const=0;

  virtual std::vector<std::array<int,2> > getEdges();
  virtual float getVol(const std::vector<Vector3f> & X);
private:
  ///@brief nodal indices
  std::vector<int> n;

};
#endif