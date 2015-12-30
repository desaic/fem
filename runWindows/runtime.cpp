#include "runtime.hpp"

int main(int argc, char* argv[])
{
  const char * filename = "config.txt";
  ConfigFile conf;
  conf.load(filename);
  int dim = conf.getInt("dim");
  if (dim == 2){
    run2D(conf);
  }
  else if (dim == 3){
    run3D(conf);
  }
  return 0;
}
