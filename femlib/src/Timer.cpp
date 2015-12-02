#include "Timer.hpp"

#include<iostream>

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

float Timer::getSeconds()
{
  return (t1 - t0) / (float)(CLOCKS_PER_SEC);
}

float Timer::getMilliseconds()
{
  return (t1 - t0) / (CLOCKS_PER_SEC/1000.0f);
}

clock_t Timer::getClocks()
{
  return t1 - t0;
}

void Timer::startWall()
{
  clock_gettime(CLOCK_MONOTONIC, &start_t);
}

void Timer::endWall()
{
  clock_gettime(CLOCK_MONOTONIC, &finish_t);
}

double Timer::getSecondsWall()
{
  double elapsed = finish_t.tv_sec - start_t.tv_sec;
  elapsed += (finish_t.tv_nsec - start_t.tv_nsec)/1.0e9;
  return elapsed;
}
