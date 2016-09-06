#include "MDS_SMACOF.h"
#include <cmath>

double mean(const Eigen::MatrixXd & v) // mean of elements
{
  double x = 0;
  for (int i = 0; i<v.rows(); i++)
  {
    for (int j = 0; j<v.cols(); j++)
    {
      x += v(i, j);
    }
  }
  return x / v.rows() / v.cols();
}


void EuclideanDistanceMatrix(const Eigen::MatrixXd & X, Eigen::MatrixXd & D)
{
  int i, j, k;
  double temp;
  
  if (X.rows() != D.rows() || X.rows() != D.cols())
  {
    printf("Invalid distance matrix dimension.\n");
    exit(1);
  }

  for (i = 0; i < D.rows(); i++) D(i, i) = 0.0;

  for (i = 0; i<D.rows() - 1; i++)
  {
    for (j = i + 1; j<D.cols(); j++)
    {
      temp = 0;
      for (k = 0; k<X.cols(); k++)
      {
        temp += std::pow(X(i, k) - X(j, k), 2);
      }
      D(i, j) = std::sqrt(temp);
    }
  }

  for (i = 1; i<D.rows(); i++)
  {
    for (j = 0; j<i; j++)
    {
      D(i, j) = D(j, i);
    }
  }
}

Eigen::MatrixXd MDS_SMACOF(const Eigen::MatrixXd & D, const Eigen::MatrixXd & X0, int dim, int iter)
{
  if (D.rows() != D.cols())
  {
    printf("Input distance matrix to MDS is not square.\n");
    exit(1);
  }
  if (dim<1)
  {
    printf("Invalid dimension for MDS.\n");
    exit(1);
  }
  if (iter<1)
  {
    printf("Invalid number of iterations for MDS.\n");
    exit(1);
  }

  Eigen::MatrixXd X;

  // with initialization
  if (X0.rows() == D.rows())
  {
    if (X0.cols() != dim)
    {
      printf("Input initialization to MDS has invalid dimension.\n");
      return X;
    }
    X = X0;
  }
  // without initialization
  else
  {
    X = Eigen::MatrixXd::Random(D.rows(), dim).array() - 0.5;
    double D_mean = mean(D); // mean value of distance matrix
    X = (0.1*D_mean / (1.0 / 3.0*sqrt((double)dim))) * X; // before this step, mean distance is 1/3*sqrt(d)
  }


  Eigen::MatrixXd Z = X;
  Eigen::MatrixXd D_ = Eigen::MatrixXd::Zero(D.rows(), D.cols());
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(D.rows(), D.cols());
  int i, j, k;
  double temp;
  double EPSILON = 1e-6;
  EuclideanDistanceMatrix(X, D_);

  printf("MDS iteration:");
  for (int it = 0; it<iter; it++) // iterations
  {
    if (it % 10 == 0) printf("\n");
    printf("%3d  ", it + 1);

    // B = calc_B(D_,D);
    for (i = 0; i<D.rows(); i++)
    {
      for (j = 0; j<D.cols(); j++)
      {
        if (i == j || fabs(D_(i, j))<EPSILON)
        {
          B(i, j)=0.0;
        }
        else
        {
          B(i, j) = -D(i, j) / D_(i, j);
        }
      }
    }

    for (j = 0; j<D.cols(); j++)
    {
      temp = 0;
      for (i = 0; i<D.rows(); i++)
      {
        temp += B(i, j);
      }
      B(j, j) =-temp;
    }

    // X = B*Z/size(D,1);
    for (i = 0; i<X.rows(); i++)
    {
      for (j = 0; j<X.cols(); j++)
      {
        temp = 0;
        for (k = 0; k<B.cols(); k++)
        {
          temp += (B(i, k)*Z(k, j));
        }
        X(i, j) =  temp / (double)D.rows();
      }
    }

    // D_ = calc_D (X);
    EuclideanDistanceMatrix(X, D_);

    Z = X;
  }

  printf("\n");
  return X;
}
