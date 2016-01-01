#include "SamplesGeneratorImpl.h"


int main(int argc, char* argv[])
{
  int dim=3;
  SamplesGeneratorImpl samplesGenerator(dim);
  samplesGenerator.setOutputDirectory("..//..//Output_2//");
  //samplesGenerator.setOutputDirectory("..//..//Output//PoissonRatio_Subdiv8_GrowAllFirst//");

  int status = samplesGenerator.run();
  return status;
  return 0;
}







