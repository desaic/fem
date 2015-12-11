#include "assert.h"
#include "cfgUtilities.h"
#include "ConfigFile.hpp"
#include "FileUtil.hpp"

#include <cmath>
#include <iostream>

int main(int argc, char * argv[]){
  const char * filename = "config.txt";
  if(argc>1){
    filename = argv[1];
  }
  ConfigFile conf;
  conf.load(filename);

  return 0;
}

void binaryToText(ConfigFile & conf)
{
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
  std::string matAssignmentFile = "../data/level16_SMC_matAssignments.bin";//fileRootName + "_matAssignments" + fileExtension;

  success = cfgUtil::readBinary<int>(matAssignmentFile, materialAssignments);
  if(!success){
    std::cout<<"Can't read "<< matAssignmentFile<<"\n";
  }

  FileUtilOut out("microstructure.txt");

  int len = (int)materialAssignments[0].size();
  std::cout<< len <<"\n";
  out.out<<materialAssignments.size()<<" "<<len<<"\n";
  for(unsigned int ii = 0;ii<materialAssignments.size(); ii++){
    for(int jj = 0; jj<len; jj++){
      out.out<<materialAssignments[ii][jj]<<" ";
    }
    out.out<<"\n";
  }

  //load and save tensors in plain text
  //  std::string tensorFile = fileRootName + "_elasticityTensors" + fileExtension;
  //  success = cfgUtil::readBinary<float>( tensorFile, tensors);
  //  if(!success){
  //    std::cout<<"Can't read "<< tensorFile<<"\n";
  //  }
  //  FileUtilOut tensorOut("elasticity.txt");
  //  int tensorLen = (int)tensors[0].size();
  //  tensorOut.out<<tensors.size()<<" "<<tensorLen<<"\n";
  //  for(unsigned int ii = 0;ii<materialAssignments.size(); ii++){
  //    for(int jj = 0; jj<tensorLen; jj++){
  //      std::cout<<tensors[ii][jj]<<" ";
  //      tensorOut.out<<tensors[ii][jj]<<" ";
  //    }
  //    tensorOut.out<<"\n";
  //  }
}
