#include "runtime.hpp"
#include "cfgUtilities.h"

///@brief convert binary to text for plotting in matlabs
void binaryParamToText(const ConfigFile & conf);

void loadMicrostructuresBin(const ConfigFile & conf, std::vector<std::vector<int> > & st);

int main(int argc, char* argv[])
{
  const char * filename = "config.txt";
  ConfigFile conf;
  conf.load(filename);

  if (conf.hasOpt("paramfile")){
    binaryParamToText(conf);
  }

  std::vector<std::vector<int> > st;
  if (conf.hasOpt("structurebin")){
    loadMicrostructuresBin(conf, st);
  }

  int dim = conf.getInt("dim");
  if (dim == 2){
    run2D(conf);
  }
  else if (dim == 3){
    run3D(conf);
  }
  return 0;
}

void binaryParamToText(const ConfigFile & conf)
{
  std::string filename = conf.dir + "/" + conf.getString("paramfile");
  bool success;
  std::vector<float> params;
  success = cfgUtil::readBinary<float>(filename, params);
  if (!success){
    std::cout << "Can't read " << filename << "\n";
  }
  std::cout << "nparams: " << params.size() << "\n";
  std::string outname("out.txt");
  if (conf.hasOpt("paramfileout")){
    outname = conf.getString("paramfileout");
  }
  std::ofstream out(outname);
  out << params.size() << "\n";
  for (unsigned int ii = 0; ii < params.size(); ii++){
    out << params[ii] << "\n";
  }
}

void loadMicrostructuresBin(const ConfigFile & conf, std::vector<std::vector<int> > & st)
{
  std::string filename = conf.dir + "/" + conf.getString("structurebin");
  bool success;
  success = cfgUtil::readBinary<int>(filename, st);
  if (!success){
    std::cout << "Can't read " << filename << "\n";
  }
  std::cout << "#structures: " << st.size() << "\n";
}
