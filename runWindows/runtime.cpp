#include "runtime.hpp"
#include "Render.hpp"
#include "World.hpp"
#include "EigenUtil.hpp"
#include "Element.hpp"
#include "cfgUtilities.h"
#include "Homogenize3D.hpp"
#include "ElementRegGrid.hpp"

///@brief convert binary to text for plotting in matlabs
void binaryParamToText(const ConfigFile & conf);

void loadMicrostructuresBin(const ConfigFile & conf, std::vector<std::vector<int> > & st);

void run3Dh(ConfigFile & conf)
{
  //testEigenUtil();
  Homogenize3D h;
  int N = 4;
  ElementRegGrid * em = new ElementRegGrid(N,N,N);
  h.gridSize = std::vector<int>(3, N);
  h.em = em;
  h.init();
  h.solve();
  for (size_t i = 0; i < em->e.size(); i++){
    em->e[i]->color = h.distribution[i] * Eigen::Vector3f(1, 1, 1);
  }
  bool render = true;
  if (render){
    Render render;
    World * world = new World();
    world->em.push_back(em);
    world->u = &h.u;
    //world->fe = &h->externalForce;
    render.init(world);
    render.loop();
  }
}

int main(int argc, char* argv[])
{
  const char * filename = "config.txt";
  ConfigFile conf;
  conf.load(filename);

  if (conf.hasOpt("paramfile")){
    binaryParamToText(conf);
  }

  int dim = conf.getInt("dim");
  if (dim == 2){
    run2D(conf);
  }
  else if (dim == 3){
    //run3D(conf);
    run3Dh(conf);
  }
  return 0;
}

void binaryParamToText(const ConfigFile & conf)
{
  std::string filename = conf.dir + "/" + conf.getString("paramfile");
  bool success;
  std::vector<float> paramOut;
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
  //float m = 0;
  //int mi=0;
  for (unsigned int ii = 0; ii < params.size(); ii++){
    //if (ii % 5 == 0){
    //  paramOut.push_back(params[ii]);
    //  paramOut.push_back(params[ii + 1]);
    //  paramOut.push_back(params[ii + 3]);
    //  paramOut.push_back(params[ii + 4]);
    //}
    //if (ii % 4 == 0){
    //  if (params[ii] > m){
    //    m = params[ii];
    //    mi = ii;
    //  }
    //}
    out << params[ii] << "\n";
  }
  //std::cout << m << " " << mi << "\n";
  //cfgUtil::writeBinary<float>("cubic-2d.bin", paramOut);
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
