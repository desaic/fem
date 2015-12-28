#include "runtime.hpp"

int main(int argc, char* argv[])
{
  const char * filename = "config.txt";
  ConfigFile conf;
  conf.load(filename);  
  run2D(conf);
  return 0;
}
