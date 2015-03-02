#include "ConfigFile.hpp"
#include "Render.hpp"
#include "World.hpp"
#include "Element.hpp"
#include "ElementHier.hpp"
#include "ElementRegGrid.hpp"
#include "ElementMeshHier.hpp"
#include "StepperGrad.hpp"
#include "StepperNewton.hpp"
#include "NewtonCuda.hpp"
#include "LinSolveCusp.hpp"
#include "AdmmCPU.hpp"
#include "AdmmNoSpring.hpp"
#include "NewtonCuda.hpp"
//#include "IpoptStepper.hpp"
#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"
#include "StrainLin.hpp"
#include "StrainCorotLin.hpp"
#include "UnitTests.hpp"

#include <iostream>
#include <iomanip>
void runTest()
{
  //ElementCoarseTest();
  stiffnessTest(2);
  //testCoarseDefGrad();
	//forceTest(0);
  testCG();
  //cudaLinTest();
  system("pause");
  exit(0);
}

void testHierForce(const ConfigFile & conf)
{
  int nlevel = 2;
  int nx = 4, ny = 4, nz = 4;
  ElementRegGrid * grid = new ElementRegGrid(nx, ny, nz);
  std::vector<StrainLin> ene(2);
  ene[0].param[0] = 10000;
  ene[0].param[1] = 100000;
  ene[1].param[0] = 100000;
  ene[1].param[1] = 1000000;
  std::vector<MaterialQuad> material(ene.size());
  for (unsigned int ii = 0; ii < material.size(); ii++){
    for (unsigned int jj = 0; jj < material[ii].e.size(); jj++){
      material[ii].e[jj] = &ene[ii];
    }
    grid->m.push_back(&material[ii]);
  }

  //assign some materials
  for (int ii = 0; ii < nx; ii++){
    for (int jj = 0; jj < ny; jj++){
      for (int kk = 0; kk < nz; kk++){
        int eidx = grid->GetEleInd(ii, jj, kk);
        if (jj < ny / 2){
          //  em->me[eidx] = 1;
        }
      }
    }
  }
  Vector3f ff(1, -4, 0);
  for (int ii = 0; ii<nx; ii++){
    for (int jj = 0; jj<nz; jj++){
      int eidx = grid->GetEleInd(ii, 0, jj);
      int aa[4] = { 0, 1, 4, 5 };
      for (int kk = 0; kk<4; kk++){
        int vidx = grid->e[eidx]->at(aa[kk]);
        grid->fixed[vidx] = 1;
      }

      eidx = grid->GetEleInd(ii, ny - 1, jj);
      int bb[4] = { 2, 3, 6, 7 };
      for (int kk = 0; kk<4; kk++){
        int vidx = grid->e[eidx]->at(bb[kk]);
        grid->fe[vidx] += ff;
      }
    }
  }

  ElementMeshHier * em = new ElementMeshHier();
  em->m.push_back(grid);
  em->buildHier(nlevel);
  ElementMesh * cm = em->m[nlevel];
  for (unsigned int ii = 0; ii < em->m[1]->e.size(); ii++){
    ElementHier * fine = (ElementHier*)(em->m[1]->e[ii]);
    std::cout << fine->parent << "\n";
    for (unsigned int jj = 0; jj < fine->Xn.size(); jj++){
      std::cout << fine->Xn[jj][0] << " " << fine->Xn[jj][1] << " " << fine->Xn[jj][2] << " | ";
    }
    std::cout << "\n";
  }

  for (unsigned int ii = 0; ii < cm->X.size(); ii++){
    std::cout << cm->X[ii][0] << " " << cm->X[ii][1] << " " << cm->X[ii][2] << " | ";
    std::cout << cm->fe[ii][0] << " " << cm->fe[ii][1] << " " << cm->fe[ii][2] << " | ";
    std::cout << cm->fixed[ii] << "\n";
  }

  Vector3f offset(0.1f, 0.2f, 0.0f);
  for (unsigned int level = 0; level < em->m.size(); level++){
    float dx = em->m[level]->X[1][2] - em->m[level]->X[0][2];
    for (unsigned int vi = 0; vi < em->m[level]->x.size(); vi++){
      if (em->m[level]->fixed[vi]){
        continue;
      }
      em->m[level]->x[vi] += dx*offset * em->m[level]->x[vi][1];
    }
    em->updateDefGrad(level);
  }

  //sanity check
  //for (unsigned int level = 0; level < em->m.size(); level++){
  //  float E = 0;
  //  for (unsigned int ei = 0; ei < em->m[level]->e.size(); ei++){
  //    E += em->getEnergy(level, ei);
  //  }
  //  std::cout << "E: " << E << "\n";
  //}

  //numerical differencing force.
  float h = 0.001f;
  for (unsigned int level = 0; level < em->m.size(); level++){
    float E = 0;
    std::cout << "====== " << level << " ======\n";
    std::vector<Vector3f> f = em->getForce(level);
    //numerical diff force
    std::vector<Vector3f> fnum(f.size());

    std::cout << std::scientific;
    std::cout << std::setprecision(3);

    ElementMesh * m = em->m[level];

    for (unsigned int vi = 0; vi < m->x.size(); vi++){
      for (int dim = 0; dim < 3; dim++){
        m->x[vi][dim] += h;
        em->updateDefGrad(level);
        float Ep, Em;
        Ep = em->getEnergy();
        m->x[vi][dim] -= 2 * h;
        em->updateDefGrad(level);
        Em = em->getEnergy();
        m->x[vi][dim] += h;
        fnum[vi][dim] = (Ep - Em) / (2 * h);
      }
    }
    float maxDiff = 0;
    for (unsigned int ii = 0; ii < f.size(); ii++){
      std::cout << f[ii][0] << " " << f[ii][1] << " " << f[ii][2] << " | ";
      std::cout << fnum[ii][0] << " " << fnum[ii][1] << " " << fnum[ii][2] << "\n";
      Vector3f diff = f[ii] - fnum[ii];
      float dd = diff.abs();
      if (dd > maxDiff){
        maxDiff = dd;
      }
    }
    std::cout << "max diff norm: " << maxDiff << "\n";
    std::cout << "E: " << E << "\n";
  }
  system("pause");
}

void testHierStiffness(const ConfigFile & conf)
{
  int nlevel = 2;
  int nx = 4, ny = 4, nz = 4;
  ElementRegGrid * grid = new ElementRegGrid(nx, ny, nz);
  std::vector<StrainLin> ene(2);
  ene[0].param[0] = 10000;
  ene[0].param[1] = 100000;
  ene[1].param[0] = 100000;
  ene[1].param[1] = 1000000;
  std::vector<MaterialQuad> material(ene.size());
  for (unsigned int ii = 0; ii < material.size(); ii++){
    for (unsigned int jj = 0; jj < material[ii].e.size(); jj++){
      material[ii].e[jj] = &ene[ii];
    }
    grid->m.push_back(&material[ii]);
  }

  //assign some materials
  for (int ii = 0; ii < nx; ii++){
    for (int jj = 0; jj < ny; jj++){
      for (int kk = 0; kk < nz; kk++){
        int eidx = grid->GetEleInd(ii, jj, kk);
        if (jj < ny / 2){
          //  em->me[eidx] = 1;
        }
      }
    }
  }

  Vector3f ff(1, -4, 0);
  for (int ii = 0; ii<nx; ii++){
    for (int jj = 0; jj<nz; jj++){
      int eidx = grid->GetEleInd(ii, 0, jj);
      int aa[4] = { 0, 1, 4, 5 };
      for (int kk = 0; kk<4; kk++){
        int vidx = grid->e[eidx]->at(aa[kk]);
        grid->fixed[vidx] = 1;
      }

      eidx = grid->GetEleInd(ii, ny - 1, jj);
      int bb[4] = { 2, 3, 6, 7 };
      for (int kk = 0; kk<4; kk++){
        int vidx = grid->e[eidx]->at(bb[kk]);
        grid->fe[vidx] += ff;
      }
    }
  }

  ElementMeshHier * em = new ElementMeshHier();
  em->m.push_back(grid);
  em->buildHier(nlevel);
  ElementMesh * cm = em->m[nlevel];

  //apply some displacement
  Vector3f offset(0.1f, 0.2f, 0.0f);
  for (unsigned int level = 0; level < em->m.size(); level++){
    float dx = em->m[level]->X[1][2] - em->m[level]->X[0][2];
    for (unsigned int vi = 0; vi < em->m[level]->x.size(); vi++){
      if (em->m[level]->fixed[vi]){
        continue;
      }
      em->m[level]->x[vi] += dx*offset * em->m[level]->x[vi][1];
    }
    em->updateDefGrad(level);
  }

  float h = 0.001;
  std::cout << "Check stiffness analytical vs numerical\n\n";
  //numerical differencing stiffness
  for (unsigned int level = 0; level < em->m.size(); level++){
    int eIdx = 0;
    float E = 0;
    std::cout << "====== " << level << " ======\n";
    MatrixXf K = em->getStiffness(level, eIdx);
    //numerical diff force
    MatrixXf Knum(K.mm, K.nn);

    std::cout << std::scientific;
    std::cout << std::setprecision(3);

    ElementMesh * m = em->m[level];
    Element * e = m->e[eIdx];
    for (unsigned int vi = 0; vi < e->nV(); vi++){
      for (int dim = 0; dim < 3; dim++){
        m->x[e->at(vi)][dim] += h;
        em->updateDefGrad(level);
        std::vector<Vector3f> fp, fm;
        fp = em->getForce(level, eIdx);
        m->x[e->at(vi)][dim] -= 2 * h;
        em->updateDefGrad(level);
        fm = em->getForce(level, eIdx);
        m->x[e->at(vi)][dim] += h;
        for (int vj = 0; vj < e->nV(); vj++){
          Vector3f df = (0.5 / h) * (fp[vj] - fm[vj]);
          for (int dimj = 0; dimj < 3; dimj++){
            Knum(3 * vi + dim, 3 * vj + dimj) = df[dimj];
          }
        }
      }
    }
    float maxDiff = 0;
    K.print(std::cout);
    std::cout << "\n----------------\n\n";
    Knum.print(std::cout);
    for (int ii = 0; ii < K.mm; ii++){
      for (int jj = 0; jj < K.nn; jj++){
        float diff = std::abs(K(ii, jj) - Knum(ii, jj));
        if (diff > maxDiff){
          maxDiff = diff;
        }
      }
    }
    std::cout << "max diff norm: " << maxDiff << "\n";
    std::cout << "E: " << E << "\n";
  }

  system("pause");
}

int runHier(const ConfigFile & conf)
{
  //int nx = 4, ny = 16, nz = 4;
  return 0;
}

int main(int argc, char* argv[])
{
  const char * filename = "config.txt";
  ConfigFile conf;
  conf.load(filename);
  if (conf.getBool("test")){
    runTest();
  }
  else if (conf.getBool("testHierForce")){
    testHierForce(conf);
    return 0;
  }
  else if (conf.getBool("testHierStiffness")){
    testHierStiffness(conf);
    return 0;
  }else if (conf.getBool("hier")){
    int ret = runHier(conf);
    return ret;
  }
  
  int refine = 1;
  if (conf.hasOpt("refine")){
    refine = conf.getInt("refine");
  }
  int res = (int)std::pow(2, refine);
  int nx = res, ny=4*res, nz=res;
  //int nx = 32, ny = 80, nz = 32;
  //int nx = 16, ny = 40, nz = 16;
  Vector3f ff(5, -10, 0);
  //per element pushing force
  ff = (1.0f / (nx*nz)) * ff;

  ElementRegGrid * em = new ElementRegGrid(nx,ny,nz);
  std::vector<StrainEneNeo> ene(2);
  //std::vector<StrainCorotLin> ene(2);
  //std::vector<StrainLin> ene(2);

  ene[0].param[0] = 10000;
  ene[0].param[1] = 100000;
  ene[1].param[0] = 100000;
  ene[1].param[1] = 1000000;
  std::vector<MaterialQuad> material(ene.size());
  for (unsigned int ii = 0; ii < material.size(); ii++){
    for (unsigned int jj = 0; jj < material[ii].e.size(); jj++){
      material[ii].e[jj] = &ene[ii];
    }
    em->addMaterial(&material[ii]);
  }

  //assign some materials
  for (int ii = 0; ii < nx; ii++){
    for (int jj = 0; jj < ny; jj++){
      for (int kk = 0; kk < nz; kk++){
        int eidx = em->GetEleInd(ii, jj, kk);
        if (jj < ny / 2){
        //  em->me[eidx] = 1;
        }
      }
    }
  }

  for(int ii = 0;ii<nx;ii++){
    for(int jj =0;jj<nz;jj++){
      int eidx= em->GetEleInd(ii,0,jj);
      int aa[4] = {0,1,4,5};
      for(int kk = 0;kk<4;kk++){
        int vidx =em->e[eidx]->at(aa[kk]);
        em->fixed[vidx] = 1;
      }

      eidx= em->GetEleInd(ii,ny-1,jj);
      int bb[4] = {2,3,6,7};
      for(int kk = 0;kk<4;kk++){
        int vidx = em->e[eidx]->at(bb[kk]);
        em->fe[vidx] += ff;
      }
    }
  }
  em->check();

  World * world = new World();
  bool renderOpt = true;
  if (conf.hasOpt("render") && conf.getBool("render") == false){
    renderOpt = false;
  }
  if (renderOpt){
    world->em.push_back(em);
  }
  Stepper * stepper = 0;
  int nSteps = 10000;
  if (conf.hasOpt("nSteps")){
    nSteps = conf.getInt("nSteps");
  }
  std::string stepperType = conf.getString("stepper");
  if (stepperType == "newton"){
    stepper = new StepperNewton();
  }
  else if (stepperType == "newtonCuda"){
    stepper = new NewtonCuda();
  }
  else if (stepperType == "ipopt"){
    //not yet implemented
    //  stepper = new IpoptStepper();
  }
  else if (stepperType == "admmCPU"){
    stepper = new AdmmCPU();
    ((AdmmCPU*)stepper)->ro0 = 100;
  }
  else if (stepperType == "grad"){
    stepper = new StepperGrad();
  }
  else{
    stepper = new AdmmNoSpring();
  }
  //stepper->rmRigid = true;
  
  stepper->nSteps = nSteps;
  stepper->init(em);
  stepper->launchThread();
  Render render;
  world->stepper = stepper;
  render.init(world);
  render.loop();
	return 0;
}
