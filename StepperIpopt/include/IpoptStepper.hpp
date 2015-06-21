#ifndef IPOPTSTEPPER_HPP
#define IPOPTSTEPPER_HPP

#include "Stepper.hpp"
#include "IpoptStepper.hpp"
#include "IpoptInterface.hpp"
#include "ElementMesh.hpp"
#include "IpIpoptApplication.hpp"

class ElementMesh;

class IpoptStepper:public Stepper
{
public:
  int NSteps;  
  virtual int oneStep();
  virtual void init(ElementMesh * _m);
  IpoptStepper();
  virtual ~IpoptStepper();

  Ipopt::SmartPtr<IpoptInterface> problem;

  Ipopt::SmartPtr<Ipopt::IpoptApplication> app;

};
#endif
