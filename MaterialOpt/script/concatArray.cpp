#include "cfgUtilities.h"
#include <iostream>
int main(int argc, char * argv[]){
  int nFile = 

  std::string fileRootName("../data/Level16");
  std::string fileExtension(".bin");
  if(conf.hasOpt("fileRootName")){
    fileRootName = conf.getString("fileRootName");
  }
  if(conf.hasOpt("fileExtension")){
    fileExtension = conf.getString("fileExtension");
  }
  std::vector<std::vector<int> > materialAssignments;
  std::vector<std::vector<float> > tensors;
  bool success ;
  std::string matAssignmentFile = fileRootName + "_matAssignments" + fileExtension;
  std::string tensorFile = fileRootName + "_elasticityTensors" + fileExtension;
  success = cfgUtil::readBinary<int>(matAssignmentFile, materialAssignments);

  return 0;
}
