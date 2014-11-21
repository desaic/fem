#ifndef IPOPTSTEPPER_HPP
#define IPOPTSTEPPER_HPP

#include "Stepper.hpp"

class ElementMesh;

class IpoptStepper:public Stepper
{
public:
  int NSteps;  
  void step(ElementMesh * mesh);
  IpoptStepper();
  virtual ~IpoptStepper();
};
#endif
