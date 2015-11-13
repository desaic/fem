#include "SamplesGeneratorImpl.h"


int main(int argc, char* argv[])
{
  int dim=2;
  SamplesGeneratorImpl samplesGenerator(dim);
  //samplesGenerator.setOutputDirectory("..//..//Output//");
  samplesGenerator.setOutputDirectory("..//..//Output//PoissonRatio_Subdiv8_GrowAllFirst//");

  int status = samplesGenerator.run();
  return status;
  return 0;
}







