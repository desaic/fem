#include "NewtonCuda.hpp"
#include "ElementMesh.hpp"
#include "femError.hpp"
#include "ArrayUtil.hpp"
#include "Timer.hpp"
#include <iostream>
NewtonCuda::NewtonCuda() :dx_tol(1e-5f), h(1.0f), linIter(10000)
{
  solver.nIter = linIter;
  m_fixRigidMotion = false;
}

void NewtonCuda::init(ElementMesh * _m)
{
  m = _m;
  std::vector<int> I, J;
  _m->stiffnessPattern(I, J, false, m_fixRigidMotion);
  solver.init(I, J);
}

int NewtonCuda::oneStep()
{

  //Timer timer;

  std::vector<Vector3f> force = m->getForce();
  float E = m->getEnergy();

  std::vector<float> val;
  //timer.start();
  if (m_fixRigidMotion)
  {
    m->getStiffnessSparse(val, false,true,true);
  }
  else
  {
    m->getStiffnessSparse(val, false,true);
  }
  //timer.end();
  //std::cout << "assembly time: " << timer.getSeconds() << "\n";

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
  if (m_fixRigidMotion)
  {
    for(int ii =0;ii<6;ii++){
      bb.push_back(0.f);
    }
  }
 
  //timer.start();
  solver.solve(val, &(bb[0]));
  //timer.end();
  //std::cout << "lin solve time: " << timer.getSeconds() << "\n";

  //timer.start();
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
  //std::cout << "total disp: " << totalMag << " ,h: " << h << "\n";
  std::cout<< "E: " << E << "\n";
  //line search
  std::vector<Vector3f> x0 = m->x;
  h = 1.0f;
  float E1 = E;
  bool valid = true;
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
        valid = false;
        break;
      }
    }
    else{
      break;
    }
  }
  
  //timer.end();
  //std::cout << "line search time: " << timer.getSeconds() << "\n";
  if (!valid){
    return 1;
  }
  return 0;
}
