#ifndef ELEMENT2D_HPP
#define ELEMENT2D_HPP

#include <array>
#include <vector>

#include "vecmath.h"
#include "cfgDefs.h"

class Element2D
{
public:
  
  Element2D(int _n=0);
  Element2D(const Element2D & e);
  Element2D(const std::vector<int> &inodes);

  ///@brief number of vertices
  virtual int nV() const;

  //@brief number of faces. Default returns 0 (not meaningful).
  virtual int nF() const;


  int & operator[](int idx);
  int operator[](int idx)const;
  int at(int idx)const;
  void resize(int size);

  ///@param X Reference vertex positions.
  ///@param x global array of Deformed vertex positions.
  Matrix2S defGrad(Vector2S p,
    const std::vector<Vector2S> & X, const std::vector<Vector2S> & x) const;

  ///@param p natural coordinates.
  ///@return interpolation weight for each dof.
  virtual std::vector<cfgScalar> shapeFun(const Vector2S & p) const=0;

  ///@param ii index of basis function.
  ///@param xx point in natural coordinates.
  ///@param X global array of rest positions.
  virtual Vector2S shapeFunGrad(int ii, const Vector2S & xx,
    const std::vector<Vector2S> & X) const=0;

  virtual std::vector<std::array<int,2> > getEdges();
  
  virtual cfgScalar getVol(const std::vector<Vector2S> & X);

  virtual Vector2S getDisp(const Vector2S & p, const std::vector<Vector2S> & X,
    const std::vector<Vector2S>x);

  const std::vector<int> & getNodeIndices()const{ return n; }

  virtual MatrixXS getMatrixB(const Vector2S & xx, const std::vector<Vector2S> & X)=0;

  //for rendering
  Vector3S color;
private:

  ///@brief nodal indices
  std::vector<int> n;

};
#endif