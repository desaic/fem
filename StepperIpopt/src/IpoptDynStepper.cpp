#include "IpoptDynStepper.hpp"

#include "Eigen/Dense"

#include <time.h>
#include <fstream>
#include <iostream>
#include <Eigen/Sparse>

using namespace Ipopt;
int IpoptDynStepper::oneStep()
{
  // Ask Ipopt to solve the problem

  ApplicationReturnStatus status = app->OptimizeTNLP(problem);

  if (status == Solve_Succeeded) {
    std::cout << std::endl << std::endl << ": *** The problem solved!" << std::endl;
    return 0;
  }
  else {
    std::cout << std::endl << std::endl << ":*** The problem FAILED!" << std::endl;
    return -1;
  }
}

void IpoptDynStepper::init(ElementMesh * _m)
{
  std::cout<<"Init ipopt dyn solver\n";
  m = _m;
  problem = new DynProblem();
  m->computeMass();
  app = IpoptApplicationFactory();
  app->Options()->SetNumericValue("tol", 0.5);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
//  app->Options()->SetStringValue("hessian_approximation","limited-memory");
//  app->Options()->SetStringValue("derivative_test","second-order");
//  app->Options()->SetNumericValue("derivative_test_perturbation",0.01);
  app->Options()->SetStringValue("hessian_constant","yes");
  app->Options()->SetIntegerValue("max_iter", NSteps);
//steps each mesh separately

 //   char outfile[32];
//    snprintf(outfile,32,"out%lu",ii);
 //   app->Options()->SetStringValue("output_file", outfile);
    problem->ele = m;
    ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Solve_Succeeded) {
      std::cout << std::endl << std::endl << "*** Error during initialization!" << std::endl;
      return;
    }

    std::cout<<"end init ipoptstepper\n";
}

IpoptDynStepper::IpoptDynStepper():NSteps(400)
{
}

IpoptDynStepper::~IpoptDynStepper()
{
}

bool DynProblem::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag,
                            IndexStyleEnum& index_style)
{
  std::cout<<"init \n";
  //3 dimensions for each vertex
  n = 3*ele->x.size();
  //no constraints
  //will implement later
  m = 0;
  nnz_jac_g = 0;

  std::vector<int> I,J;
  ele->stiffnessPattern(I,J,true);
  // the Hessian's number of nonzeros
  //This matrix is symmetric - specify the lower diagonal only
  nnz_h_lag = J.size();
  std::cout<<"nnz "<<nnz_h_lag<<"\n";
  //zero based
  index_style = TNLP::C_STYLE;

  //h*h*K
  float h = ele->dt;
  //get full stiffness matrix
  A = ele->getStiffnessSparse();
  A=h*h*A;
//  std::cout<<A<<"\n";
  //add lumped mass
  //column first
  for(int ii = 0; ii<A.cols(); ii++){
    int vidx = ii/3;
    for (Eigen::SparseMatrix<float>::InnerIterator it(A, ii); it; ++it){
//      std::cout<<ii<<" "<<it.row()<<"\n";
      if(ii==it.row()){
//        std::cout<<it.value()<<"\n";
        it.valueRef() += ele->mass[vidx];
//        std::cout<<it.value()<<"\n";
      }
    }
  }
  b.resize(3*ele->x.size());
  std::vector<Vector3f> forces = ele->getForce();
  for(unsigned int ii =0; ii<forces.size(); ii++){
    //add momentum and gravity scaled by h
    Vector3f f = (h*h)*(forces[ii] + ele->mass[ii] * ele->G);
    f += h*ele->mass[ii] * ele->v[ii];
    for(int jj =0 ; jj<3; jj++){
      b[3*ii+jj] = f[jj];
    }
  }
  return true;
}

bool
DynProblem::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
  assert(n==(int)3*ele->x.size());
  assert(m==0);
  for(size_t ii = 0;ii<ele->x.size();ii++){
    //x and z component have lower bound -10 (arbitrary)
    //y component has lower bound 0
    x_l[3*ii] = -10;
    x_l[3*ii+1] = 0;//-10;
    x_l[3*ii+2] = -10;

    //all variables upper bound 10 (arbitrary)
    x_u[3*ii] = 10;
    x_u[3*ii+1] = 10;
    x_u[3*ii+2] = 10;

    if(ele->fixed[ii]){
      x_l[3*ii]   = ele->x[ii][0];
      x_l[3*ii+1] = ele->x[ii][1];
      x_l[3*ii+2] = ele->x[ii][2];

      //all variables upper bound 10 (arbitrary)
      x_u[3*ii]   = ele->x[ii][0];
      x_u[3*ii+1] = ele->x[ii][1];
      x_u[3*ii+2] = ele->x[ii][2];
    }
  }

  //did not implement other constraints

  return true;
}

bool
DynProblem::get_starting_point(Index n, bool init_x, Number* x,
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
//      std::cout<<x[3*ii+jj]<<"\n";
    }
  }
  return true;
}

bool
DynProblem::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  Eigen::VectorXf delta(3*ele->x.size());
  for(size_t ii = 0 ;ii<ele->x.size();ii++){
    for(int jj = 0;jj<3;jj++){
      delta[3*ii+jj] = x[3*ii+jj] - ele->x[ii][jj];
//      std::cout<<x[3*ii+jj]<<" ";
    }
//    std::cout<<"\n";
  }

  obj_value = delta.dot(0.5*A*delta - b);
  return true;
}

bool
DynProblem::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  Eigen::VectorXf delta(3*ele->x.size());
  for(size_t ii = 0 ;ii<ele->x.size();ii++){
    for(int jj = 0;jj<3;jj++){
      delta[3*ii+jj] = x[3*ii+jj] - ele->x[ii][jj];
    }
  }

  Eigen::VectorXf grad = A*delta - b;
  for(Index ii = 0; ii<n; ii++){
    grad_f[ii] = grad[ii];
  }
  return true;
}

// return the value of the constraints: g(x)
bool
DynProblem::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  //0 at the moment

  return true;
}

bool
DynProblem::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
  if (values == NULL) {
  }
  else {
  }
  return true;
}

//return the structure or values of the hessian.
//This matrix is symmetric - specify the lower diagonal only.
//not used if LBFGS.
bool
DynProblem::eval_h(Index n, const Number* x, bool new_x,
                       Number obj_factor, Index m, const Number* lambda,
                       bool new_lambda, Index nele_hess, Index* iRow,
                       Index* jCol, Number* values)
{
  if(values==NULL){
    std::vector<int> I,J;
    //use only trianglar part
    ele->stiffnessPattern(I,J,true);
    int Jidx = 0;
    for(unsigned int ii = 1;ii<I.size();ii++){
      int len = I[ii] - I[ii-1];
      for(int jj =0 ; jj<len; jj++){
        iRow[Jidx] = ii-1;
        jCol[Jidx] = J[Jidx];
//        std::cout<<iRow[Jidx]<<" "<<jCol[Jidx]<<"\n";
        Jidx++;
      }
    }
    std::cout<<Jidx<<" nnz pattern\n";
  }else{

//    std::vector<int> I,J;
    int k =0 ;
    for(int ii = 0; ii<A.cols(); ii++){
      for (Eigen::SparseMatrix<float>::InnerIterator it(A, ii); it; ++it){
        //only copy triangular part
        if(it.col()>it.row()){
          continue;
        }
        values[k] = obj_factor*it.value();
        k++;
      }
    }

  }
  return true;
}


void
DynProblem::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
          const IpoptData* ip_data,
          IpoptCalculatedQuantities* ip_cq)
{
  // here is where we would store the solution to variables, or write to a file, etc
  // so we could use the solution.
  for (size_t ii=0; ii<ele->x.size(); ii++) {
    for(int jj = 0;jj<3;jj++){
      ele->v[ii][jj] = (1.0/ele->dt) * (x[3*ii+jj] - ele->x[ii][jj]);
      ele->x[ii][jj] = x[3*ii+jj];
    }
  }

  std::cout<<"Step: "<<frameCnt<<"\n";
  frameCnt++;
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

DynProblem::DynProblem():
  ele(0),frameCnt(0)
{
}

DynProblem::~DynProblem()
{}

