#include "UnitTests.hpp"
#include <stdlib.h>
void runTest()
{
  //ElementCoarseTest();
  stiffnessTest(0);
	forceTest(0);
//  testCG();
  //cudaLinTest();
  system("pause");
  exit(0);
}

int main(int argc, char* argv[]){
    runTest();
}
