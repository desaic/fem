#include "cfgUtilities.h"
#include <iostream>
#include <string>

int main(int argc, char * argv[]){
  int nFile = 11;

  std::string fileRootName("../../data/Level16_elasticityTensors/level16_");
  std::string fileExtension("_elasticityTensors.bin");
  
  std::vector<std::vector<float> > all;
  std::vector<std::vector<float> > tensors;
  bool success ;
  for(int ii = 0; ii<nFile; ii++){
    std::string tensorFile = fileRootName +std::to_string(ii) + fileExtension;
    success = cfgUtil::readBinary<float>(tensorFile, tensors);
    all.insert(all.end(), tensors.begin(), tensors.end());
  }
  
  std::string outfile("level16_elasticityTensors.bin");
  cfgUtil::writeBinary<float>(outfile, all);
  
  return 0;
}
