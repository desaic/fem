#include "SamplesGeneratorImpl.h"


int main(int argc, char* argv[])
{
  int dim=2;
  SamplesGeneratorImpl samplesGenerator(dim);
  samplesGenerator.setOutputDirectory("..//..//Output//");

  int status = samplesGenerator.run();
  return status;
  return 0;
}







