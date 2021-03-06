#include "ElementHex.hpp"
#include <iostream>

using namespace Eigen;

const int cubeEdges[12][2] =
{
  {0,1},
  {0,2},
  {1,3},
  {2,3},

  {4,5},
  {4,6},
  {5,7},
  {6,7},
  
  {0,4},
  {1,5},
  {2,6},
  {3,7}
};

static const int sw[8][3] =
{{-1,-1,-1},
 {-1,-1, 1}, 
 {-1, 1,-1},
 {-1, 1, 1},
 { 1,-1,-1},
 { 1,-1, 1},
 { 1, 1,-1},
 { 1, 1, 1}
};

const int faces[6][4]={
		{ 0, 1, 3, 2 },
		{ 4, 6, 7, 5 },
		{ 0, 4, 5, 1 },
		{ 2, 3, 7, 6 },
		{ 0, 2, 6, 4 },
		{ 1, 5, 7, 3 }
};

///@brief face normals
const int facen[6][3]={
	{-1, 0, 0},
	{ 1, 0, 0},
	{ 0,-1, 0},
	{ 0, 1, 0},
	{ 0, 0,-1},
	{ 0, 0, 1}
};

Vector3S ElementHex::natCoord(const Vector3S & p, const std::vector<Vector3S> & X)
{
  Vector3S n = p - X[at(0)];
  cfgScalar size = (X[at(7)][0] - X[at(0)][0]);
  n = (2.0f / size) * n - Vector3S(1,1,1);
  return n;
}

std::vector<cfgScalar>
ElementHex::shapeFun(const Vector3S & p)const
{
  std::vector<cfgScalar> weights(8);
  for(int ii = 0;ii<nV();ii++){
    weights[ii] = (1.0f/8) * (1+sw[ii][0]*p[0])
      *(1+sw[ii][1]*p[1]) *(1+sw[ii][2]*p[2]) ;
  }
  return weights;
}

Vector3S
ElementHex::shapeFunGrad(int ii, const Vector3S & xx,
                                 const std::vector<Vector3S> & X) const
{
  Vector3S size=4*(X[at(7)] - X[at(0)]);
//  std::cout<<size[0]<<"\n";

  Vector3S grad;
  size[0] = 1.0f/(size[0]);
  size[1] = 1.0f/(size[1]);
  size[2] = 1.0f/(size[2]);
  grad[0] = sw[ii][0] * size[0] * (1 + sw[ii][1] * xx[1]) * (1 + sw[ii][2] * xx[2]);
  grad[1] = sw[ii][1] * size[1] * (1 + sw[ii][0] * xx[0]) * (1 + sw[ii][2] * xx[2]);
  grad[2] = sw[ii][2] * size[2] * (1 + sw[ii][0] * xx[0]) * (1 + sw[ii][1] * xx[1]);

  return grad;
}

std::vector<std::array<int,2> >
ElementHex::getEdges()
{
  int nEdge = 12;
  std::vector<std::array<int,2> >  edges(nEdge);
  for(int ii=0; ii<nEdge; ii++){
    edges[ii][0] = cubeEdges[ii][0];
    edges[ii][1] = cubeEdges[ii][1];
  }
  return edges;
}

Vector3S ElementHex::facePt(int fi, const Vector3S & x)
{
	Vector3S p;
	switch (fi){
	case 0:
		p[0] = -1;
		p[1] = x[0];
		p[2] = x[1];
		break;
	case 1:
		p[0] = 1;
		p[1] = x[0];
		p[2] = x[1];
		break;
	case 2:
		p[1] = -1;
		p[0] = x[0];
		p[2] = x[1];
		break;
	case 3:
		p[1] = 1;
		p[0] = x[0];
		p[2] = x[1];
		break;
	case 4:
		p[2] = -1;
		p[0] = x[0];
		p[1] = x[1];
		break;
	case 5:
		p[2] = 1;
		p[0] = x[0];
		p[1] = x[1];
		break;
	default:
		break;
	}
	return p;
}

MatrixXS ElementHex::NMatrix(int fi)
{
	MatrixXS N= MatrixXS::Zero(3,6);
	for (int ii = 0; ii < 3; ii++){
		N(ii, ii) = (cfgScalar)facen[fi][ii];
	}
	N(0, 3) = (cfgScalar)facen[fi][1];
	N(0, 5) = (cfgScalar)facen[fi][2];

	N(1, 3) = (cfgScalar)facen[fi][0];
	N(1, 4) = (cfgScalar)facen[fi][2];

	N(2, 4) = (cfgScalar)facen[fi][1];
	N(2, 5) = (cfgScalar)facen[fi][0];
	return N;
}

MatrixXS ElementHex::HMatrix(const Vector3S & xx)
{
	std::vector<cfgScalar> w = shapeFun(xx);
	MatrixXS H = MatrixXS::Zero(3,3*nV());
	for (int ii = 0; ii < nV(); ii++){
		for (int jj = 0; jj < 3; jj++){
			H(jj, 3 * ii + jj) = w[ii];
		}
	}
	return H;
}

MatrixXS ElementHex::BMatrix(const Vector3S & xx, const std::vector<Vector3S>X)
{
	MatrixXS B = MatrixXS::Zero(6, 3 * nV());
	for (int ii = 0; ii < nV(); ii++){
		int col = 3 * ii;
		Vector3S dN = shapeFunGrad(ii, xx, X);
		B(0, col  ) = dN[0];
		B(1, col+1) = dN[1];
		B(2, col+2) = dN[2];

		B(3, col    ) = dN[1];
		B(3, col + 1) = dN[0];

		B(4, col + 1) = dN[2];
		B(4, col + 2) = dN[1];

		B(5, col    ) = dN[2];
		B(5, col + 2) = dN[0];
	}
	return B;
}

cfgScalar ElementHex::getVol(const std::vector<Vector3S> & X)
{
  Vector3S size = X[at(7)] - X[at(0)];
  cfgScalar vol = size[0] * size[1] * size[2];
//  std::cout<<vol<<"\n";
  return vol;
}

ElementHex::ElementHex():Element(8)
{}

ElementHex::ElementHex(const ElementHex & e) :Element(e)
{}

ElementHex::ElementHex(const std::vector<int> &iNodes)
  :Element(iNodes)
{
}
