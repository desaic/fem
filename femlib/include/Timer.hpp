#ifndef TIMER_HPP
#define TIMER_HPP
#include <ctime>

class Timer
{
public:
  Timer();
  void start();
  ///@return -1 if overflow happened
  int end();
  float getSeconds();
  float getMilliseconds();
  clock_t getClocks();

  void   startWall();
  void   endWall();
  double getSecondsWall();

private:
  clock_t t0, t1;

  struct timespec start_t;
  struct timespec finish_t;

};

#endif
