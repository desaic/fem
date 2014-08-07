#ifndef ADMMCPU_HPP
#define ADMMCPU_HPP
#include "Stepper.hpp"
class AdmmCPU:public Stepper
{
public:
  void step(ElementMesh * m);
  AdmmCPU();
  ///@brief a copy of vertex positions for each element.
  ///Vertices are ordered to the individual element
  std::vector<std::vector<Vector3f> > u;

  ///@lagrange multiplier for each vertex
  std::vector<std::vector<Vector3f> > y;
  
  ///@brief magic
  std::vector<float> ro;
  ///@brief number of elements sharing vertex i.
  std::vector<int> N;
  ///@brief concensus variables
  std::vector<Vector3f> Z;
  float tol;
  ///@brief tolerance for x variables.
  float xtol;
  void minimizeElement(ElementMesh * mesh, Element * ele, int eIdx);

  float getEnergy(ElementMesh *eMesh, Element * ele, int eIdx);
  
  std::vector<Vector3f>getForces(ElementMesh * mesh, Element * ele,
                         int eIdx);
  
  MatrixXd stiffness(ElementMesh *mesh, Element * ele, int eIdx);
    
private:
  void initVar(ElementMesh * e);
};

#endif