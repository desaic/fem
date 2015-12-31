#ifndef RUNTIME_HPP
#define RUNTIME_HPP

#include "ConfigFile.hpp"
#include "FEM2DFun.hpp"
#include "FEM3DFun.hpp"

#include <fstream>

void run2D(const ConfigFile & conf);
void run3D(const ConfigFile & conf);

///@brief sparse matrix vector product.
Eigen::VectorXd
spmv(const std::vector<int> & I, const std::vector<int> & J, const std::vector<double> & val,
const std::vector<double> x, bool triangle = true);

///@brief verify consistency of f and df using central differencing
void check_df(RealFun * fun, const Eigen::VectorXd & x0, double h);

void loadArr2D(std::vector<std::vector< int> > &arr, std::string filename);
void loadText(const ConfigFile & conf, std::vector<std::vector<double> > & materialAssignments);

///@brief verify that linear static simulation is working.
///Checks K*u - f.
///Also checks df/dx = K.
///Numerical differencing is very slow. It computes all N*N entries in K, where N = Number of elements.
///@TODO: Move code to UnitTests.
void check_sim(FEM2DFun * fem, const Eigen::VectorXd & x);

///@brief finds a step size the decreases value of fun.
///@return <0 if cannot find a positive step size.
///@param x0 starting point.
///@param dir search direction.
///@param h appropriate step size if there is one. Always positive.
int lineSearch(RealFun * fun, Eigen::VectorXd & x0, const Eigen::VectorXd & dir, double & h);

///@param nSteps maximum number of steps.
void gradientDescent(RealFun * fun, Eigen::VectorXd & x0, int nSteps);

void loadIntBinary(const ConfigFile & conf, std::vector<std::vector<double> > & materialAssignments);

///@brief compute material properties
void computeMat(FEM2DFun * fem, const ConfigFile & conf);

///@brief optimize material distribution for each starting point.
void optMat(FEM2DFun * fem, int nSteps);

void optMat3D(FEM3DFun * fem, int nSteps);

#endif