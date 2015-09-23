#include "Quadrature2D.h"

const cfgScalar Gauss2Pt[2]=
{
 (cfgScalar)-1/std::sqrt((cfgScalar)3),
  (cfgScalar)1/std::sqrt((cfgScalar)3)
};

Quadrature2D makeGauss2();
Quadrature2D makeUniform4();

const Quadrature2D Quadrature2D::Gauss2=makeGauss2();
const Quadrature2D Quadrature2D::Uniform4=makeUniform4();

Quadrature2D makeGauss2()
{
  Quadrature2D q;
  q.x.resize(4);
  q.w.resize(4,0.25f);
  q.x[0] = Vector2S(Gauss2Pt[0],Gauss2Pt[0]);
  q.x[1] = Vector2S(Gauss2Pt[0],Gauss2Pt[1]);
  q.x[2] = Vector2S(Gauss2Pt[1],Gauss2Pt[0]);
  q.x[3] = Vector2S(Gauss2Pt[1],Gauss2Pt[1]);
  return q;
}

Quadrature2D makeUniform4()
{
  Quadrature2D q;
  q.w.resize(16, (cfgScalar)1/(cfgScalar)16);
  int nPt = 4;
  cfgScalar dx = (cfgScalar)2/nPt;
  cfgScalar x0 = -1;
  for(int ii = 0; ii<nPt; ii++){
    for(int jj = 0; jj<nPt; jj++){
      Vector2S p(x0+dx*(ii+(cfgScalar)0.5),x0+dx*(jj+(cfgScalar)0.5));
      q.x.push_back(p); 
    }
  }
  return q;
}

Quadrature2D::Quadrature2D()
{}

