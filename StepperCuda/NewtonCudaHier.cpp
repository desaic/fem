#include "NewtonCudaHier.hpp"
#include "NewtonCuda.hpp"
#include "ElementMesh.hpp"
#include "ElementMeshHier.hpp"
#include "femError.hpp"
#include "ArrayUtil.hpp"
#include "Timer.hpp"
#include <iostream>
NewtonCudaHier::NewtonCudaHier() :dx_tol(1e-5f), h(1.0f), linIter(10000)
{
  solver.nIter = linIter;
}

void NewtonCudaHier::init(ElementMesh * _m)
{
  m = _m;
  ElementMeshHier * emh = (ElementMeshHier*)m;
  level = (int)(emh->m.size() - 1);

  std::vector<int> I, J;
  emh->getStiffnessPattern(level, I, J, false);
  solver.init(I, J);
}

int NewtonCudaHier::oneStep()
{
  ElementMeshHier * emh = (ElementMeshHier *)(m);
  ElementMesh * em = emh->m[level];
  //Timer timer;
  std::vector<Vector3f> force = emh->getForce(level);
  float E = emh->getEnergy(level);

  std::vector<float> val;
  //timer.start();
  emh->getStiffness(level, val, false, true);
  //timer.end();
  //std::cout << "assembly time: " << timer.getSeconds() << "\n";

  int ndof = 3 * (int)em->x.size();
  std::vector<float> bb(ndof);

  for (unsigned int ii = 0; ii<em->x.size(); ii++){
    for (int jj = 0; jj<3; jj++){
      int row = 3 * ii + jj;
      if (em->fixed[ii]){
        bb[row] = 0;
      }
      else{
        bb[row] = force[ii][jj];
      }
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
  for (unsigned int ii = 0; ii<em->x.size(); ii++){
    for (int jj = 0; jj<3; jj++){
      force[ii][jj] = (float)bb[3 * ii + jj];
    }
  }
  //std::cout << "total disp: " << totalMag << " ,h: " << h << "\n";
  std::cout << "E: " << E << "\n";
  //line search
  std::vector<Vector3f> x0 = em->x;
  h = 1.0f;
  float E1 = E;
  bool valid = true;
  while (1){
    if (totalMag * h<dx_tol){
      valid = false;
      break;
    }
    em->x = x0;
    addmul(em->x, h, force);
    E1 = emh->getEnergy(level);
    if (E1>=E || fem_error){
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
    if (level > 0){
      level--;
      solver.dealloc();
      std::vector<int> I, J;
      emh->getStiffnessPattern(level, I, J, false);
      solver.init(I, J);
      return 0;
    }
    else{
      return -1;
    }
  }
  return 0;
}
