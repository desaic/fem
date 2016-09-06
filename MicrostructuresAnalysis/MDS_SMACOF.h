/// \file MDS_SMACOF.h
#include <Eigen/Dense>

double mean(const Eigen::MatrixXd & v); // mean of elements

void EuclideanDistanceMatrix(const Eigen::MatrixXd & X, Eigen::MatrixXd & D);

/// \brief http://tosca.cs.technion.ac.il/
/// \param X0 guess initial embedding parameters. 
/// Pass a matrix with wrong number of rows, e.g. 1x1, if no guess is provided.
Eigen::MatrixXd MDS_SMACOF(const Eigen::MatrixXd & D, const Eigen::MatrixXd & X0, int dim, int iter);
