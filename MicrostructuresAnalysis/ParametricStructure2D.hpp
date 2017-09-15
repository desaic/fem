#ifndef PARAMETRIC_STRUCTURE_2D_HPP
#define PARAMETRIC_STRUCTURE_2D_HPP

#include "MicrostructureParam2D.hpp"
#include <vector>

//microstructure described by a set of parameters.
class ParametricStructure2D{
public:

  ParametricStructure2D(){ 
    s.resize(2, 0);
  }

  ParametricStructure2D(int Nx, int Ny);

  void resize(int Nx, int Ny);

  //shape parameters.
  std::vector<double> param;

  //updated 2D structure after calling eval.
  std::vector<double> s;

  //2D grid size of s.
  std::vector<int> gridSize;
  virtual void eval(const std::vector<double> & p){ 
    param = p;
  }
};

//negative poisson's ratio parametric structure
class ParametricNPR :public ParametricStructure2D{
public:
  ParametricNPR();
  //square in center
  Rectangle2D rec_center;
  //branch off the center square
  Rectangle2D rec_branch;
  //short beams
  Rectangle2D rec_beam;
  //cross like structures
  Rectangle2D rec_cross;
  Ellipse2D rec_cross_tip;
  //distribute parameters into geometric primitives.
  void distributeParam();
  
  //read variable parameters from geometric primitives
  void extractParam();

  /// \brief meaning of 7 parameters:
  /// center rectangle size with equal x and y size.
  /// branch length
  /// branch width
  /// beam length
  /// beam width
  /// cross length
  /// cross width
  void eval(const std::vector<double> & p);

};

void fitTarget(ParametricStructure2D * st , const std::vector<double> & target,
  const std::vector<int> & gridSize);

void test2DParam();

#endif