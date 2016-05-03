#include "MicrostructureOpt.h"
#include "ElasticHexFEM.h"
#include "Homogenization.h"

static const int d = 3;
Type_Define_VectorD_Outside_Class(d); Type_Define_VectorDi_Outside_Class(d); Type_Define_MatrixD_Outside_Class(d);

struct MicrostructureOpt
{
  Lattice<d> lattice;
  ElasticHexFEM<d> fem;
  Homogenization<d> * hom;
};

void saveDisplacements(const MicrostructureOpt & opt);

void initOpt(MicrostructureOpt & o, PARSE_ARGS & parse_args)
{
  const int s = parse_args.Get_Integer_Value("-s");
  const int solver_type = parse_args.Get_Integer_Value("-solver");
  const int processor_type = parse_args.Get_Integer_Value("-proc");
  const bool write_residual = parse_args.Get_Option_Value("-wr");
  const int num_mg_level = parse_args.Get_Integer_Value("-mg_lv");
  T length = (T)1; VectorDi cell_counts = VectorDi::Ones()*(s-1);
  o.lattice = Lattice<d>(cell_counts, length / (T)(cell_counts[0]+1));
  o.fem.Initialize(o.lattice);
  o.fem.solver_type = ElasticHexFEM<d>::SolverType(solver_type);
  o.fem.processor_type = ElasticHexFEM<d>::ProcessorType(processor_type);
  o.fem.des_iter_num = parse_args.Get_Integer_Value("-n_d");
  o.fem.asc_iter_num = parse_args.Get_Integer_Value("-n_a");
  o.fem.v_cycle_iter_num = parse_args.Get_Integer_Value("-n_c");
  o.fem.krylov_iter_num = parse_args.Get_Integer_Value("-krylov_iter_num");
  int material_pattern = parse_args.Get_Integer_Value("-mt");
  o.fem.use_implicit_psi_P = true;
  o.fem.num_multigrid_level = num_mg_level;

  ////reinitialize materials
  T E = 1;
  T nu = 0.45;
  o.fem.materials.clear();
  o.fem.materials.push_back(ElasticMaterial(E, nu));
  for (LatticeIterator<d> iter(o.lattice, LatticeIterator<d>::CELL); iter.Valid(); iter.Next()){
    o.fem.material_id(iter.Coord()) = 0;
  }

  o.hom = new Homogenization<d>(o.fem);
  o.hom->Initialize();

}

void saveDisplacements(const MicrostructureOpt & opt)
{
  for (size_t i = 0; i < opt.hom->unit_test_u.size(); i++){
    std::string filename = "0m" + std::to_string(i) + ".txt";
    std::ofstream out(filename);
    out << opt.hom->unit_test_u[i].size() / d << "\n";
    for (size_t j = 0; j < opt.fem.node_num; j++){
      for (int k = 0; k < d; k++){
        T disp = opt.hom->unit_test_u[i][d*j + k];
        out << disp << " ";
      }
      out << "\n";
    }
  }
}

void computeCTensor(const MicrostructureOpt & opt)
{
  int n = 6;
  Eigen::MatrixXd C = Eigen::MatrixXd::Zero(n, n);
  VectorD domainSize = opt.lattice.domain_max - opt.lattice.domain_min; 
  T vol = domainSize.prod();
  int i = 0;
  int j = 0;
  
  for (int i = 0; i < n; i++)for (int j = i; j < n; j++){
    for (LatticeIterator<d> iter(opt.lattice, LatticeIterator<d>::CELL); iter.Valid(); iter.Next()){
      VectorDi cell = iter.Coord(); int mat_id = opt.fem.material_id(cell);
      if (mat_id != -1){
        VectorX cell_u_i; opt.fem.Compute_Cell_Displacement_Helper(opt.hom->unit_test_u[i], cell, cell_u_i);
        VectorX cell_u_j; opt.fem.Compute_Cell_Displacement_Helper(opt.hom->unit_test_u[j], cell, cell_u_j);

        if (i == 0){
          for (int k = 0; k < cell_u_j.size(); k++){
            std::cout << k << " " << cell_u_j[k] << "\n";
          }
          std::cout << "============\n";
        }

        T e = cell_u_i.dot(opt.fem.Ke[mat_id] * cell_u_j);
        C(i, j) += e;
      }
    }
  }
  int input;
  T a = C(0, 1) / C(0, 0); T poisson = a / ((T)1 + a);
  T youngs=C(0,1)*((T)1+poisson)*((T)1-(T)2*poisson)/poisson;
  std::cout << youngs<<" "<<poisson << "\n";
  std::cin >> input;
}

void TestMicrostructureOpt(PARSE_ARGS & parse_args)
{
  std::string conffile("../config.txt");
  ConfigFile conf;
  conf.load(conffile.c_str());
  MicrostructureOpt o;

  //initOpt(o, parse_args);
  //o.hom->Homogenize();
  //computeCTensor(o);
  run3D(conf, parse_args);
}