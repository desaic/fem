#include "cfgUtilities.h"
#include "ConfigFile.hpp"
#include "EigenUtil.hpp"
#include "FEM3DFun.hpp"
#include "FileUtil.hpp"
#include "PARSE_ARGS.h"
#include "PiecewiseConstant3D.hpp"
#include "PiecewiseConstantSym3D.hpp"
#include "PiecewiseConstantCubic3D.hpp"

#include <algorithm>
#include <iomanip>
#include <thread>

extern std::ofstream logfile;
static const int d = 3;
Type_Define_VectorD_Outside_Class(d); Type_Define_VectorDi_Outside_Class(d); Type_Define_MatrixD_Outside_Class(d);

void loadIntBinary(const ConfigFile & conf, std::vector<std::vector<double> > & materialAssignments);

struct Opt3DArgs
{
  FEM3DFun * fem;
  const ConfigFile * conf;
};

void shrinkVector(Eigen::VectorXd & x, float shrink_ratio)
{
  for (int ii = 0; ii < x.size(); ii++){
    x[ii] = 0.5 + shrink_ratio * (x[ii] - 0.5);
  }
}

void optMat3D(Opt3DArgs * arg)
{
  int nSteps = 500;
  const ConfigFile * conf = arg->conf;
  FEM3DFun * fem = arg->fem;

  RealField * field = fem->field;
  
    //for (int ii = 0; ii < x1.size(); ii++){
  //  x1[ii] = 0.5 + (rand() / (float)RAND_MAX - 0.5) * 0.2;
  //}
 
  std::vector<std::vector<double> > structures;
  //std::vector<int> idx;
  //cfgUtil::readBinary<int>("../data/boundaryPoints.bin", idx);
  loadIntBinary(*conf, structures);
  std::ofstream matStruct("struct.txt");
  double shrinkRatio = 0.3;
  for (unsigned int si = 0; si < structures.size(); si++){
    std::cout << si << "\n";
    std::vector<double> s3d = structures[si];
    Eigen::VectorXd x1 = fem->param;
    //copy 8x8x8 corner to 3D structure.
    std::vector<int> paramSize = fem->gridSize;
    for (int dim = 0; dim < (int)paramSize.size(); dim++){
      paramSize[dim] /= 2;
    }
    for (int ix = 0; ix < paramSize[0]; ix++){
      for (int iy = 0; iy < paramSize[1]; iy++){
        for (int iz = 0; iz < paramSize[2]; iz++){
          int linearIdx = ix * paramSize[1] * paramSize[2] + iy *paramSize[2] + iz;
          //int inputIdx = ix / 2 * paramSize[1] * paramSize[2] + iy / 2 * paramSize[2] + iz / 2;
          int inputIdx = ix * fem->gridSize[1] * fem->gridSize[2] + iy * fem->gridSize[2] + iz;
          x1[linearIdx] = 1;// s3d[inputIdx];
        }
      }
    }
    clampVector(x1, fem->lowerBounds, fem->upperBounds);
    fem->setParam(x1);
    T a = fem->G(0, 1) / fem->G(0, 0); T poisson = a / ((T)1 + a);
    T youngs = fem->G(0, 1)*((T)1 + poisson)*((T)1 - (T)2 * poisson) / poisson;
    std::cout << fem->G(0, 0) << " " << fem->G(1, 0) << " " << fem->G(2, 0) << " " << fem->G(3, 3) << "\n";
    std::cout << "rho: " << fem->density << ", E: " << youngs << ", nu: " << poisson << "\n";
    system("pause");
    //double val = fem->f();
    //fem->m0 = 0.5*fem->density;
    ////scale mass term to roughly displacement term.
    //fem->mw = 0.1 * fem->G(0, 0) / fem->density;
    //fem->G0 = fem->G;
    ////poisson's ratio objective
    //fem->G0(1, 0) = 0.5 * fem->G0(0, 0);
    //fem->G0(2, 0) = 0.5 * fem->G0(0, 0);

    ////std::cout << fem->G << "\n";
    ////for test only look at the first displacement.
    //for (int ii = 1; ii < fem->wG.size(); ii++){
    //  fem->wG(ii) = 0;
    //}
    //shrinkVector(x1, shrinkRatio);
    ////logfile.open("log3d.txt");
    ////check_df(fem, x1, 1e-3);
    //gradientDescent(fem, x1, nSteps);
    //for (unsigned int ii = 0; ii < fem->distribution.size(); ii++){
    //  matStruct << fem->distribution[ii] << " ";
    //}
    //matStruct << "\n";
  }
  system("pause");
}

void initFEM(FEM3DFun & fem, PARSE_ARGS & parse_args , int gridRes)
{
  const int solver_type = parse_args.Get_Integer_Value("-solver");
  const int processor_type = parse_args.Get_Integer_Value("-proc");
  const bool write_residual = parse_args.Get_Option_Value("-wr");
  const int num_mg_level = parse_args.Get_Integer_Value("-mg_lv");
  //initialize without one layer ghost cell.
  //size edge length to 1 including ghost cell
  T length = (T)1; VectorDi cell_counts = VectorDi::Ones()*(gridRes - 1);
  fem.lattice = Lattice<d>(cell_counts, length / (T)(cell_counts[0] + 1));
  fem.fem.Initialize(fem.lattice);
  fem.fem.solver_type = ElasticHexFEM<d>::SolverType(solver_type);
  fem.fem.processor_type = ElasticHexFEM<d>::ProcessorType(processor_type);
  fem.fem.des_iter_num = parse_args.Get_Integer_Value("-n_d");
  fem.fem.asc_iter_num = parse_args.Get_Integer_Value("-n_a");
  fem.fem.v_cycle_iter_num = parse_args.Get_Integer_Value("-n_c");
  fem.fem.krylov_iter_num = parse_args.Get_Integer_Value("-krylov_iter_num");
  int material_pattern = parse_args.Get_Integer_Value("-mt");
  fem.fem.use_implicit_psi_P = true;
  fem.fem.num_multigrid_level = num_mg_level;

  ////reinitialize materials
  T E = 1;
  T nu = 0.45;
  fem.fem.materials.clear();
  fem.fem.materials.push_back(ElasticMaterial(E, nu));
  for (LatticeIterator<d> iter(fem.lattice, LatticeIterator<d>::CELL); iter.Valid(); iter.Next()){
    fem.fem.material_id(iter.Coord()) = 0;
  }

  fem.fem.u.setZero(); fem.fem.f.setZero(); fem.fem.psi_D.clear(); fem.fem.nodal_forces.clear();
  fem.fem.Update_Matrices();
  fem.Ke0 = fem.fem.Ke;
  int n = 6;
  fem.u.resize(6);
  int nVert = std::pow(gridRes+1, d);
  for (int i = 0; i < n; i++){
    fem.u[i].resize(d * nVert);
  }

  Eigen::Vector3i varSize( gridRes,  gridRes, gridRes);
  fem.hom = new Homogenization<d>(fem.fem);
  fem.hom->Initialize();
  fem.fem.variable_coef = new Field<T, d>();
  fem.fem.variable_coef->Resize(varSize); fem.fem.variable_coef->Fill((T)0);

}

void run3D(const ConfigFile & conf, PARSE_ARGS & parse_args)
{
  bool render = conf.getBool("render");
  std::string task = conf.getString("task");

  int gridres = 32;
  double maxStep = 1.0;
  if (conf.hasOpt("gridres")){
    gridres = conf.getInt("gridres");
  }
  if (conf.hasOpt("maxStep")){
    maxStep = conf.getFloat("maxStep");
  }
  int nx = gridres;
  int ny = gridres;
  int nz = gridres;
  int nSteps =500;
  Eigen::VectorXd x0;

  FEM3DFun * fem = new FEM3DFun();
  initFEM(*fem, parse_args, gridres);
  //PiecewiseConstant3D * field = new PiecewiseConstant3D();
  //PiecewiseConstantSym3D * field = new PiecewiseConstantSym3D();
  PiecewiseConstantCubic3D * field = new PiecewiseConstantCubic3D();
  field->allocate(nx / 2, ny / 2, nz / 2 );
  //field->allocate(nx , ny , nz);
  fem->lowerBounds = 1e-3 * Eigen::VectorXd::Ones(field->param.size());
  fem->upperBounds = Eigen::VectorXd::Ones(field->param.size());
  fem->param = field->param;
  fem->distribution.resize(nx*ny*nz);
  float fx = 1;

  fem->field = field;
  fem->gridSize[0] = nx;
  fem->gridSize[1] = ny;
  fem->gridSize[2] = nz;
  fem->m0 = 0.4;
  fem->mw = 0.1;
  field->param = Eigen::VectorXd::Ones(field->param.size());
  Opt3DArgs * arg = new Opt3DArgs;
  arg->fem = fem;
  arg->conf = &conf;
  optMat3D(arg);
}


void loadIntBinary(const ConfigFile & conf, std::vector<std::vector<double> > & materialAssignments)
{
  std::string fileName("../data/level16_matAssignments.bin");
  if (conf.hasOpt("structurebin")){
    fileName = conf.getString("structurebin");
  }
  bool success;
  std::string matAssignmentFile = conf.dir + "/" + fileName;
  std::vector<std::vector<int> > intArr;
  success = cfgUtil::readBinary<int>(matAssignmentFile, intArr);
  if (!success){
    std::cout << "Can't read " << matAssignmentFile << "\n";
    return;
  }
  materialAssignments.resize(intArr.size());
  for (unsigned int ii = 0; ii < intArr.size(); ii++){
    materialAssignments[ii].resize(intArr[ii].size());
    for (unsigned int jj = 0; jj < materialAssignments[ii].size(); jj++){
      materialAssignments[ii][jj] = intArr[ii][jj];
    }
    //int nx = (int)(std::sqrt(materialAssignments[ii].size()) + 0.4);
    //materialAssignments[ii] = upsampleVector(materialAssignments[ii], nx, nx);
  }

}
