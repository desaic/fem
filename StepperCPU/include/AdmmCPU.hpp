#ifndef ADMMCPU_HPP
#define ADMMCPU_HPP
#include "Stepper.hpp"
#include "vecmath.h"
#include <vector>
#include <fstream>
#include "MatrixX.hpp"
class Element;
class AdmmCPU:public Stepper
{
public:
  int oneStep();
  void init(ElementMesh * _m);
  AdmmCPU();

  ~AdmmCPU();
  
  ///@brief a copy of vertex positions for each element.
  ///Vertices are ordered to the individual element
  std::vector<std::vector<Vector3f> > u;

  ///@lagrange multiplier for each vertex
  std::vector<std::vector<Vector3f> > y;
  
  ///@brief magic
  std::vector<float> ro;
  
  ///@brief concensus variables
  std::vector<Vector3f> Z;
  
  ///@temporary array for lin solve
  float * bb ;

  ///@brief maximum distance allowed between x and z
  float maxDist;

  float ro0;

  float roMult;

  float tol;
  
  ///@brief tolerance for x variables.
  float xtol;
  
  void minimizeElement(ElementMesh * mesh, Element * ele, int eIdx);

  float getEnergy(ElementMesh * eMesh, int eIdx);
  
  std::vector<Vector3f> getForces(ElementMesh * eMesh, int eIdx);
  
  MatrixXf stiffness(ElementMesh *mesh, int eIdx);

  float prevE;
  std::ofstream out;
private:
  void initVar();
};

#endif