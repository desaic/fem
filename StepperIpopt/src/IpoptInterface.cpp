#include "TimeStepper/IpoptInterface.hpp"
#include "Element/ElementMesh.hpp"
#include "Eigen/Dense"
#include "World/ConstraintVertBd.hpp"

#include <time.h>
#include <fstream>
#include <iostream>
#include <Eigen/Sparse>

using namespace Ipopt;
int frameCnt =0;
bool IpoptInterface::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag,
                            IndexStyleEnum& index_style)
{
  //3 dimensions for each vertex
  n = 3*ele->x.size();
  //no constraints
  //will implement later
  m = 0;
  nnz_jac_g = 0;

  // the Hessian's number of nonzeros
  //This matrix is symmetric - specify the lower diagonal only
  nnz_h_lag = ele->stiffness_pattern().nonZeros();
  //std::cout<<"nnz "<<nnz_h_lag<<"\n";
  //zero based
  index_style = TNLP::C_STYLE;
  return true;
}

bool
IpoptInterface::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  assert(n==(int)3*ele->x.size());
  assert(m==0);
  for(size_t ii = 0;ii<ele->x.size();ii++){
    //x and z component have lower bound -10 (arbitrary)
    //y component has lower bound 0
    x_l[3*ii] = -10;
    x_l[3*ii+1] = -10;
    x_l[3*ii+2] = -10;

    //all variables upper bound 10 (arbitrary)
    x_u[3*ii] = 10;
    x_u[3*ii+1] = 10;
    x_u[3*ii+2] = 10;
  }

  for(size_t ii = 0 ; ii<ele->vertConst.size();ii++){
    ConstraintVertBd & cc = ele->vertConst[ii];
    for(int jj = 0;jj<3;jj++){
      x_l[ 3*cc.vertIdx + jj ] = cc.l[jj];
      x_u[ 3*cc.vertIdx + jj ] = cc.u[jj];
    }
  }

  //did not implement other constraints

  return true;
}

bool
IpoptInterface::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  assert(init_x == true);
  assert(init_z == false);
  assert(init_lambda == false);

  // initialize to the given starting point
  for(size_t ii = 0;ii<ele->x.size();ii++){
    for(int jj = 0;jj<3;jj++){
      x[ii*3+jj] = ele->x[ii][jj];
    }
  }
  return true;
}

bool
IpoptInterface::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  for(size_t ii = 0 ;ii<ele->x.size();ii++){
    for(int jj = 0;jj<3;jj++){
      ele->x[ii][jj] = x[3*ii+jj];
    }
  }
  obj_value = ele->getEnergy();
  return true;
}

bool
IpoptInterface::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  for(size_t ii = 0 ;ii<ele->x.size();ii++){
    for(int jj = 0;jj<3;jj++){
      ele->x[ii][jj] = x[3*ii+jj];
    }
  }

  std::vector<Eigen::Vector3f> forces = ele->GetForces();
  assert((int)3*forces.size() == n);

  for(size_t ii = 0;ii<forces.size();ii++){
    for(int jj = 0;jj<3;jj++){
      grad_f[ii*3+jj] = -forces[ii][jj];
    }
  }
  return true;
}

// return the value of the constraints: g(x)
bool
IpoptInterface::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  //0 at the moment

  return true;
}

bool
IpoptInterface::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  if (values == NULL) {
  }
  else {
  }
  return true;
}

std::vector<float>eneSeq;
std::vector<float>timeSeq;

//return the structure or values of the hessian.
//This matrix is symmetric - specify the lower diagonal only.
//not used if LBFGS.
bool
IpoptInterface::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
  clock_t tt;
  eneSeq.push_back(ele->getEnergy());
  tt = clock();
  timeSeq.push_back(tt/(CLOCKS_PER_SEC/1000.0));
  if(values==NULL){
    int cnt = 0;
    Eigen::SparseMatrix<float> mat = ele->stiffness_pattern();
    for (int k=0; k<mat.outerSize(); ++k){
      for (Eigen::SparseMatrix<float>::InnerIterator it(mat,k); it; ++it){
        iRow[cnt] = it.row();
        jCol[cnt] = it.col();
//        std::cout<<iRow[cnt] <<" "<<jCol[cnt]<<"\n";
        cnt++;
      }
    }

    assert(cnt == nele_hess);
  }else{
    ele->saveObj(frameCnt);
    frameCnt ++ ;
    int cnt=0;
    Eigen::SparseMatrix<float> mat = ele->stiffness();
    std::cout<<"nnz "<<mat.nonZeros()<<"\n";
    for (int k=0; k<mat.outerSize(); ++k){
      for (Eigen::SparseMatrix<float>::InnerIterator it(mat,k); it; ++it){
        values[cnt] = obj_factor*it.value();
        cnt++;
      }
    }
    
    for (size_t ii=0; ii<ele->x.size(); ii++) {
    for(int jj = 0;jj<3;jj++){
      ele->x[ii][jj] = x[3*ii+jj];
    }
    }
    
  }
  return true;
}


void
IpoptInterface::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
          const IpoptData* ip_data,
          IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.

  // Copy the solution to the elements.
  // This is why we need a compiler project.
  for (size_t ii=0; ii<ele->x.size(); ii++) {
    for(int jj = 0;jj<3;jj++){
      ele->x[ii][jj] = x[3*ii+jj];
    }
  }
  ele->saveObj(frameCnt);
  frameCnt ++ ;


  clock_t tt;
  eneSeq.push_back(ele->getEnergy());
  tt = clock();
  timeSeq.push_back(tt/(CLOCKS_PER_SEC/1000.0));
  std::ofstream out("converge.txt");
  for(unsigned int ii = 1;ii<eneSeq.size();ii++){
    float relTime = timeSeq[ii]-timeSeq[0];
    out<<relTime<<" "<<eneSeq[ii]<<"\n";
  }
  out.close();

#ifdef DEBUG
  std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
  for (Index i=0; i<n; i++) {
    std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
  }
  for (Index i=0; i<n; i++) {
    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
  }

  std::cout << std::endl << std::endl << "Objective value" << std::endl;
  std::cout << "f(x*) = " << obj_value << std::endl;

  std::cout << std::endl << "Final value of the constraints:" << std::endl;
  for (Index i=0; i<m ;i++) {
    std::cout << "g(" << i << ") = " << g[i] << std::endl;
  }
#endif
}

IpoptInterface::IpoptInterface():
  ele(0)
{
}

IpoptInterface::~IpoptInterface()
{}
