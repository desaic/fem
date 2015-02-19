#include "Timer.hpp"
Timer::Timer() :t0(0), t1(0)
{}

void Timer::start()
{
  t0 = std::clock();
}

int Timer::end()
{
  t1 = std::clock();
  //naive overflow check
  if (t1 < t0){
    return -1;
  }
  return 0;
}

float Timer::getMilliseconds()
{
  return (t1 - t0) / (float)(CLOCKS_PER_SEC);
}

clock_t Timer::getClocks()
{
  return t1 - t0;
}
