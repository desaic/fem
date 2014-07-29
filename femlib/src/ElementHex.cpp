#include "ElementHex.hpp"

int  cubeEdges[12][2]=
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

int sw[8][3] =
{{-1,-1,-1},
 {-1,-1, 1}, 
 {-1, 1,-1},
 {-1, 1, 1},
 { 1,-1,-1},
 { 1,-1, 1},
 { 1, 1,-1},
 { 1, 1, 1}
};

std::vector<float>
ElementHex::ShapeFun(const Vector3f & p)const
{
  std::vector<float> weights(8);
  for(int ii = 0;ii<nV();ii++){
    weights[ii] = (1.0f/8) * (1+sw[ii][0]*p[0])
      *(1+sw[ii][1]*p[1]) *(1+sw[ii][2]*p[2]) ;
  }
  return weights;
}

Vector3f
ElementHex::ShapeFunGrad(int ii, const Vector3f & xx,
                                 const std::vector<Vector3f> & X) const
{
  Vector3f size=4*(X[7] - X[0]);
	Vector3f grad;
	size[0] = 1.0f/(size[0]);
	size[1] = 1.0f/(size[1]);
	size[2] = 1.0f/(size[2]);

  grad[0] = sw[ii][0] * size[0] * (1 + sw[ii][1] * xx[1]) * (1 + sw[ii][2] * xx[2]);
  grad[1] = sw[ii][1] * size[1] * (1 + sw[ii][0] * xx[0]) * (1 + sw[ii][2] * xx[2]);
  grad[2] = sw[ii][2] * size[2] * (1 + sw[ii][0] * xx[0]) * (1 + sw[ii][1] * xx[1]);

	return grad;
}

std::vector<std::array<int,2> >
ElementHex::GetEdges()
{
  int nEdge = 12;
  std::vector<std::array<int,2> >  edges(nEdge);
  for(int ii=0; ii<nEdge; ii++){
    edges[ii][0] = cubeEdges[ii][0];
    edges[ii][1] = cubeEdges[ii][1];
  }
  return edges;
}

ElementHex::ElementHex():Element(8)
{}