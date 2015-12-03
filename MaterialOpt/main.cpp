#include "assert.h"
#include "cfgUtilities.h"
#include "ConfigFile.hpp"
#include <cmath>
#include <iostream>
int main(int argc, char * argv[]){
  const char * filename = "config.txt";
  if(argc>1){
    filename = argv[1];
  }
  ConfigFile conf;
  conf.load(filename);

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
  if(!success){
    std::cout<<"Can't read "<< matAssignmentFile<<"\n";
  }
  success = cfgUtil::readBinary<float>( tensorFile, tensors);
  if(!success){
    std::cout<<"Can't read "<< matAssignmentFile<<"\n";
  }
  for(unsigned int ii = 0;ii<materialAssignments.size(); ii++){
    int len = (int)materialAssignments[ii].size();
    std::cout<< len <<"\n";
    int nx = std::sqrt(len);
    for(int jj = 0; jj<len; jj++){
      int xi = jj/nx;
      int yi = jj%nx;
      int pxi = yi;
      int pyi = nx-xi-1;
      std::cout<<materialAssignments[ii][pxi*nx+pyi]<<" ";
      if(pxi == nx-1){
        std::cout<<"\n";
      }
    }

    len = (int)tensors[ii].size();
    for(int jj = 0; jj<len; jj++){
      std::cout<<tensors[ii][jj]<<" ";
    }
    std::cout<<"\n";
  }

  return 0;
}
