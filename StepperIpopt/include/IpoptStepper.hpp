#ifndef IPOPTSTEPPER_HPP
#define IPOPTSTEPPER_HPP
#include "TimeStepper.hpp"
class World;

class IpoptStepper:public TimeStepper
{
public:
  int NSteps;  
  void Step(World * world);
  IpoptStepper();
  virtual ~IpoptStepper();
};
#endif
