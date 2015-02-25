#include "NewtonCuda.hpp"
#include "ElementMesh.hpp"
#include "femError.hpp"
#include "ArrayUtil.hpp"

#include <iostream>
NewtonCuda::NewtonCuda() :dx_tol(1e-5f), h(1.0f), linIter(1000)
{
  solver.nIter = linIter;
}

void NewtonCuda::init(ElementMesh * _m)
{
  m = _m;
  std::vector<int> I, J;
  _m->stiffnessPattern(I, J, false);
  solver.init(I, J);
}

int NewtonCuda::oneStep()
{
  std::vector<Vector3f> force = m->getForce();
  float E = m->getEnergy();

  std::vector<float> val;
  m->getStiffnessSparse(val, false,true);

  int ndof = 3 * (int)m->x.size();
  std::vector<float> bb(ndof);

  for (unsigned int ii = 0; ii<m->x.size(); ii++){
    for (int jj = 0; jj<3; jj++){
      int row = 3 * ii + jj;
      if (m->fixed[ii]){
        bb[row] = 0;
      }
      else{
        bb[row] = force[ii][jj];
      }
    }
  }
 
  solver.solve(val, &(bb[0]));

  double totalMag = 0;
  for (int ii = 0; ii<ndof; ii++){
    totalMag += std::abs(bb[ii]);
  }
  totalMag /= ndof;
  for (unsigned int ii = 0; ii<m->x.size(); ii++){
    for (int jj = 0; jj<3; jj++){
      force[ii][jj] = (float)bb[3 * ii + jj];
    }
  }
  std::cout << "total disp: " << totalMag <<" ,h: "<<h<< ", E:"<<E<<"\n";
  //line search
  std::vector<Vector3f> x0 = m->x;
  h = 1.0f;
  float E1 = E;
  while (1){
    if (totalMag * h<dx_tol){
      return -1;
    }
    m->x = x0;
    addmul(m->x, h, force);
    E1 = m->getEnergy();
    if (E1>E || fem_error){
      fem_error = 0;
      h = 0.5f* h;
      if (h<1e-12){
        return -1;
        break;
      }
    }
    else{
      return 0;
    }
  }
  
}
