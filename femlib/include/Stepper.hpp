#ifndef STEPPER_HPP
#define STEPPER_HPP
class ElementMesh;

class Stepper{
public:
  Stepper();
  virtual void step(ElementMesh * m)=0;
  virtual ~Stepper();

  int nSteps;
};

#endif