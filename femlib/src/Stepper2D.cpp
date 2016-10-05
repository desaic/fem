#include "Stepper2D.h"
#include "Timer.hpp"

#include <atomic>
#include <iostream>

Stepper2D::Stepper2D() :nSteps(100), m(0)
{}

Stepper2D::~Stepper2D(){}

void runSim(ElementMesh2D * m, Stepper2D * Stepper2D)
{
  Timer t;
  for (int iter = 0; iter < Stepper2D->nSteps; iter++){
    std::cout << "iter: " << iter <<"\n";
    std::unique_lock<std::mutex> lck(Stepper2D->mtx);
    while (Stepper2D->state == Stepper2D::PAUSE){
      Stepper2D->cv.wait(lck);
    }
    if (Stepper2D->state == Stepper2D::SINGLE){
      Stepper2D->state = Stepper2D::PAUSE;
    }
    lck.unlock();
    
    t.start();
    int ret = Stepper2D->oneStep();
    t.end();
    float duration = t.getSeconds();
    //std::cout << "time: " << duration << "\n";
    if (ret < 0){
      break;
    }
  }
  std::cout << "Exit Stepper2D loop\n";
}

void Stepper2D::launchThread()
{
  state = PAUSE;
  thread = std::thread(runSim, m, this);
}

void Stepper2D::init(ElementMesh2D * _m)
{
  m = _m;
}

int Stepper2D::step()
{
  Timer t;
  for (int ii = 0; ii < nSteps; ii++){
    t.start();
    int ret = oneStep();
    t.end();
    float duration = t.getSeconds();
    //std::cout << ii << " time: " << duration << "\n";
    if (ret < 0){
      return ret;
    }
    else if (ret==1)
    {
      std::cout << "Failed!" << std::endl;
      return ret;
    }
  }
  return 0;
}
