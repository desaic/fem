#include "Stepper.hpp"

#include <atomic>
#include <chrono>
Stepper::Stepper() :nSteps(100), m(0)
{}

Stepper::~Stepper(){}

void runSim(ElementMesh * m, Stepper * stepper)
{
  for (int iter = 0; iter < stepper->nSteps; iter++){
    
    std::unique_lock<std::mutex> lck(stepper->mtx);
    while (stepper->state == Stepper::PAUSE){
      stepper->cv.wait(lck);
    }
    if (stepper->state == Stepper::SINGLE){
      stepper->state = Stepper::PAUSE;
    }
    lck.unlock();

    int ret = stepper->oneStep();
    if (ret < 0){
      break;
    }
  }
}

void Stepper::launchThread()
{
  state = PAUSE;
  thread = std::thread(runSim, m, this);
}

void Stepper::init(ElementMesh * _m)
{
  m = _m;
}

int Stepper::step()
{
  for (int ii = 0; ii < nSteps; ii++){
    int ret = oneStep();
    if (ret < 0){
      return ret;
    }
  }
  return 0;
}
