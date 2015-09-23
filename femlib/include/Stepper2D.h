#ifndef STEPPER_2D_HPP
#define STEPPER_2D_HPP
#include <thread>
#include <mutex>
#include <condition_variable>

class ElementMesh2D;

class Stepper2D{
public:
  enum StepperState{PAUSE, SINGLE, ALL};
  
  Stepper2D();
  virtual void init(ElementMesh2D * _m);
  virtual int step();
  virtual int oneStep() = 0;
  virtual ~Stepper2D();
  
  void launchThread();
  void notify(StepperState s);
  
  int nSteps;
  ElementMesh2D * m;

  StepperState state;
  std::thread thread;
  std::mutex mtx;
  std::condition_variable cv;
};

#endif