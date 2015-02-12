#ifndef ADMMNOSPRING_HPP
#define ADMMNOSPRING_HPP
#include <Eigen/Dense>
#include "vecmath.h"
#include "Stepper.hpp"
#include <vector>
#include <fstream>

class Element;
struct MatrixXd;
class AdmmNoSpring :public Stepper
{
public:

  int oneStep();
  void init(ElementMesh * _m);
  AdmmNoSpring();

  ~AdmmNoSpring();

  ///@brief a copy of vertex positions for each element.
  ///Vertices are ordered to the individual element
  std::vector<std::vector<Vector3f> > u;

  ///@lagrange multiplier for each vertex
  std::vector<std::vector<Vector3f> > y;
  
  ///@brief concensus variables
  std::vector<Vector3f> Z;

  ///@temporary array for lin solve
  double * bb;

  float tol;

  ///@brief tolerance for x variables.
  float xtol;

  void minimizeElement(ElementMesh * mesh, Element * ele, int eIdx);

  float getEnergy(ElementMesh * eMesh, int eIdx);

  std::vector<Vector3f> getForces(ElementMesh * eMesh, int eIdx);

  MatrixXd stiffness(ElementMesh *mesh, int eIdx);

  std::vector<Eigen::MatrixXf> T;


  float prevE;
  std::ofstream out;
};

#endif