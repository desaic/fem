#ifndef ADMM_H
#define ADMM_H
#include <vector>
#include <string>
#include "Stepper.hpp"
#include "vecmath.h"
#include "vector_types.h"
struct ADMMInfo;
class ElementMesh;
class ADMMStepper:public Stepper
{
public:
	ADMMStepper();
  ~ADMMStepper();

  void init(ElementMesh * _m);

  ///@lagrange multiplier for each copy of a vertex
  std::vector<std::vector<Vector3f> > l;
  
  int NSteps;
  int nThread;

  ///initial ro for all variables
  float ro;
  ///@brief number of elements sharing vertex i.
  std::vector<int> N;
  ///@brief concensus variables
  std::vector<Vector3f> Z;
  float tol;
  ///@brief maximum distance between x and Z in each dimension
  ///unitless. ratio with respect to element size
  float maxDist;
  ///@brief ro is increased by multiplying roMult each time
  float roMult;
  float3 * hostxx, *hostXX;
  float3 * devXX, * devxx;
  ADMMInfo * devadmm, * hostadmm;

  float  * Ehost, *Edev;
  float3 * fhost, *fdev;
  float prevE;

  int oneStep();
  void initVar();

  //used for line search
  
  //copy vec to devadmm.Z
  void setZdev(const std::vector<Vector3f> & vec);
  
  //get energy of mesh using current Z stored in devadmm, excluding ADMM terms
  float getEnergy();
  
  //get forces for the mesh excluding ADMM terms
  void getForce(std::vector<Vector3f> & ff);
  
  //check gradient and line search.
  //incredibly slow convergence.
  //not used for simulation.
  void stepGrad();

  std::string outname;
};
#endif