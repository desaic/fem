#include "AdmmCPU.hpp"
#include "Element.hpp"
#include "ElementHex.hpp"
#include "ElementMesh.hpp"
#include <time.h>
#include <fstream>

float
ADMMStepper::getEnergy(ElementMesh * eMesh, Element * ele,
                       int eIdx)
{
  float E = ele->GetEnergy(eMesh);
  for(int ii = 0;ii<ele->getNumNodes();ii++){
    Eigen::Vector3f diff = eMesh->x[ele->At(ii)]-Z[ele->At(ii)]+l[eIdx][ii];
    E += 0.5* ro * diff.dot(diff);
  }
  return E;
}

std::vector<Eigen::Vector3f>
ADMMStepper::getForces(ElementMesh * eMesh, Element * ele,
                       int eIdx)
{
  std::vector<Eigen::Vector3f> force = ele->GetNodalForces(eMesh);
  for(size_t ii = 0;ii<force.size();ii++){
    Eigen::Vector3f diff = eMesh->x[ele->At(ii)] -Z[ele->At(ii)]+l[eIdx][ii];
    //when minimizing force is negative gradient
    force[ii] -= ro * diff;
  }

  return force;
}

Eigen::MatrixXf
ADMMStepper::stiffness(ElementMesh *mesh, Element *ele, int eIdx)
{
  Eigen::MatrixXf K = ((ElementHex*)ele)->Stiffness(mesh);
  for(int ii = 0;ii<K.rows();ii++){
    K(ii,ii) += ro;
  }
  return K;
}

void ADMMStepper::testStiff(ElementMesh * mesh, Element * ele,
                            int eIdx)
{
  float E = getEnergy(mesh,ele,eIdx);
  float h = 0.0001;
  std::vector<Eigen::Vector3f>force = getForces(mesh,ele,eIdx);
  for(auto ii = 0;ii<ele->getNumNodes();ii++){
    for(int jj = 0; jj<3; jj++){
      mesh->x[ele->At(ii)][jj] -= h;
      double Eminus = getEnergy(mesh,ele,eIdx);
      mesh->x[ele->At(ii)][jj] += 2*h;
      double Eplus = getEnergy(mesh,ele,eIdx);
      std::cout<<"Energy diff: "<<Eplus-Eminus<<"\n";
      double numF = (Eplus - Eminus )/(2*h);
      std::cout<<"Analytic derivative:\n"<<-force[ii][jj];
      std::cout<<"\nCentral diff:\n"<<numF<<"\n--------\n";
    }
  }

  Eigen::MatrixXf KAna = stiffness(mesh,ele,eIdx);
  int nVar = ele->getNumNodes();
  Eigen::MatrixXf K = Eigen::MatrixXf::Zero(3*nVar,3*nVar);
  for(auto ii = 0; ii<ele->getNumNodes(); ii++){
    for(int jj = 0; jj<3; jj++){
      mesh->x[ele->At(ii)][jj] -= h;
      std::vector<Eigen::Vector3f>fminus = getForces(mesh,ele,eIdx);
      mesh->x[ele->At(ii)][jj] += 2*h;
      std::vector<Eigen::Vector3f>fplus = getForces(mesh,ele,eIdx);
      mesh->x[ele->At(ii)][jj] -=h;
      for(size_t kk = 0; kk<fminus.size(); kk++){
        fplus[kk] -= fminus[kk];
        fplus[kk] /= 2*h;
        for(int ll = 0;ll<3;ll++){
          K(3*kk+ll,3*ii+jj)=-fplus[kk][ll];
        }
      }
    }
  }
  Eigen::MatrixXf diff = K - KAna;
  std::cout<<"Num K:\n"<<KAna<<"\n";
  float maxErr = 0;
  for(int k=0; k<diff.outerSize(); ++k){
    for(auto ii = 0;ii<K.rows();ii++){
      for(auto jj = 0;jj<K.cols();jj++){
        float error = std::abs(diff(ii,jj));
        if(error>maxErr){
          maxErr = error;
        }
      }
    }
  }
  std::cout<<"max err "<<maxErr<<"\n";
}

std::vector<Eigen::Vector3f>
getEleX(ElementMesh * mesh, Element * ele)
{
  std::vector<Eigen::Vector3f> xx(ele->getNumNodes());
  for(int ii = 0;ii<ele->getNumNodes();ii++){
    xx[ii] = mesh->x[ele->At(ii)];
  }
  return xx;
}

std::vector<Eigen::Vector3f>
setEleX(ElementMesh * mesh, Element * ele, const std::vector<Eigen::Vector3f> & xx)
{
  for(int ii = 0;ii<ele->getNumNodes();ii++){
    mesh->x[ele->At(ii)] = xx[ii];
  }
}

void
ADMMStepper::minimizeElement(ElementMesh * mesh, Element * ele,
                             int eIdx)
{
 // testStiff(mesh,ele,eIdx);

  float E = 0;
  int ii = 0;
  float localh = 1;
  int NSteps = 100;

  for(ii = 0;ii<NSteps;ii++){
    E = getEnergy(mesh,ele,eIdx);
    std::vector<Eigen::Vector3f> nodalForces = getForces(mesh,ele,eIdx);
    std::vector<Eigen::Vector3f> x0(ele->getNumNodes());
    for(int ii= 0;ii<ele->getNumNodes();ii++){
      x0[ii] = mesh->x[ele->At(ii)];
    }
    Eigen::MatrixXf K = stiffness(mesh,ele,eIdx);
    Eigen::VectorXf bb = Eigen::VectorXf(K.rows());
    for(unsigned int ii = 0;ii<nodalForces.size();ii++){
      for(int jj = 0;jj<3;jj++){
        bb[3*ii+jj] = nodalForces[ii][jj];
      }
    }

    bb=K.colPivHouseholderQr().solve(bb);
//    std::cout<<"E0 "<<E<<"\n";
    nodalForces = getEleX(mesh,ele);

    float change;
    while(true){
      change = 0;
      for(int ii = 0; ii<ele->getNumNodes(); ii++){
        for(int jj = 0; jj<3; jj++){
          mesh->x[ele->At(ii)][jj] = nodalForces[ii][jj] + localh*bb(3*ii+jj);
          change += std::abs(bb[3*ii+jj]);
        }
      }
      change *= localh;
      float ene =getEnergy(mesh,ele,eIdx);
//      for(auto ii = 0;ii<mesh->x.size();ii++){
//        for(int jj = 0; jj<3; jj++){
//          std::cout<<mesh->x[ii][jj]<<" ";
//        }
//        std::cout<<"\n";
//      }
//      for(auto ii = 0;ii<mesh->x.size();ii++){
//        for(int jj = 0; jj<3; jj++){
//          std::cout<<Z[ii][jj]<<" ";
//        }
//        std::cout<<"\n";
//      }
      if(std::isnan(ene) || ene>=E){
        localh/=2;
//        std::cout<<"Step size:" <<localh<<" Energy: "<<ene<<" change: "<<change<<"\n";
        if(localh<0.00001 || change<xtol){
          return;
        }
      }else{
        break;
      }
    }
    if(change<xtol){
//      std::cout<<"Steps: "<<ii<<" "<<change<<"\n";
      break;
    }
  }
}

void ADMMStepper::initVar(ElementMesh *e)
{
  u.resize(e->elements.size());
  l.resize(u.size());
  N.assign(e->X.size(),0);
  for(size_t ii= 0; ii<e->elements.size();ii++){
    Element * ele = e->elements[ii];
    u[ii].resize(ele->getNumNodes());
    l[ii].resize(u[ii].size(),Eigen::Vector3f::Zero());
    for(int jj = 0;jj<ele->getNumNodes();jj++){
      u[ii][jj] = e->X[ele->At(jj)];
      N[ele->At(jj)] ++;
    }
  }

  for(auto ii = 0;ii< e->nForces.size();ii++){
    ForceNode * ff = e->nForces[ii];
    N[ff->index] ++;
  }

  cl.assign(e->vertConst.size(), Eigen::Vector3f::Zero());
  for(auto ii = 0; ii< e->vertConst.size(); ii++){
    N[e->vertConst[ii].vertIdx] ++;
  }

  e->x = e->X;
  Z = e->x;
}

void ADMMStepper::Step(ElementMesh * e)
{
  initVar(e);
  float E = e->getEnergy();
  float prevE = E;

  std::ofstream out("converge.txt");
  clock_t tt,tt0 ;
  tt0 = clock();

  std::vector<Eigen::Vector3f> Z_k = e->x;

  bool pause = true;
  e->u = &u;
  for(int iter = 0;iter<NSteps;iter++){
    if(pause){
      std::string str;
      std::cin>>str;
      if(str[0] == 'c'){
        pause = false;
      }
    }
    if(iter%100==0){
      pause=true;
    }

    //primal and dual residual
    float rr = 0, ss=0;
    //objective. Augmentedd lagrangian
    float obj = 0;
    //update u locally
    for(size_t ii = 0;ii<e->elements.size();ii++){
      Element * ele = e->elements[ii];
      for(size_t jj = 0;jj<ele->getNumNodes();jj++){
        e->x[ele->At(jj)] = u[ii][jj]; //Z[ele->At(jj)]
      }
      minimizeElement(e,ele, ii);

      obj += ele->GetEnergy(e);

      for(size_t jj = 0;jj<ele->getNumNodes();jj++){
        u[ii][jj] = e->x[ele->At(jj)];
      }
    }

    for(size_t ii = 0;ii<u.size();ii++){
      Element * ele = e->elements[ii];
      for(size_t jj = 0;jj<ele->getNumNodes();jj++){
        obj += ro * l[ii][jj].dot(u[ii][jj] - Z[ele->At(jj)]);
        obj += 0.5*ro*(u[ii][jj] - Z[ele->At(jj)]).squaredNorm();
      }
    }
    std::cout<<iter<<" Objective: "<<obj<<"\n";

    std::vector<Eigen::Vector3f> uf(e->nForces.size());
    //force variables solved in closed form
    for(size_t ii = 0;ii<e->nForces.size();ii++){
      ForceNode * f= e->nForces[ii];
      int vidx = f->index;
      uf[ii] = Z[vidx] + (1/ro) * f->force ;
    }

    //update z closed form
    for(size_t ii = 0;ii<e->X.size();ii++){
      Z[ii] = Eigen::Vector3f(0,0,0);
    }

    //add per element variables and multipliers
    for(size_t ii = 0;ii<u.size();ii++){
      Element * ele = e->elements[ii];
      for(size_t jj = 0;jj<ele->getNumNodes();jj++){
        Z[ele->At(jj)] += u[ii][jj]+l[ii][jj];
      }
    }

    //add force variables
    for(size_t ii = 0;ii<e->nForces.size();ii++){
      ForceNode * f= e->nForces[ii];
      int vidx = f->index;
      Z[vidx] += uf[ii];
    }

    //add constraint variables
    for(size_t ii = 0;ii<e->vertConst.size();ii++){
      Eigen::Vector3f uc = e->vertConst[ii].l;
      int vidx = e->vertConst[ii].vertIdx;
      Z[vidx] += uc + cl[ii];
    }

    //divide
    for(size_t ii = 0;ii<e->x.size();ii++){
      Z[ii] /= N[ii];
      ss += N[ii] * (Z_k[ii] - Z[ii]).squaredNorm();
    }

    //update multiplier for elements
    for(size_t ii = 0;ii<u.size();ii++){
      Element * ele = e->elements[ii];
      for(size_t jj = 0;jj<ele->getNumNodes();jj++){
        l[ii][jj] +=  u[ii][jj] - Z[ele->At(jj)];
        rr        += (u[ii][jj] - Z[ele->At(jj)]).squaredNorm();
      }
    }

    //update multiplier for constraints
    for(size_t ii = 0; ii<cl.size(); ii++){
      int vidx = e->vertConst[ii].vertIdx;
      cl[ii] += e->vertConst[ii].l - Z[vidx];
    }

    Z_k = Z;
    e->x = Z;
    E = e->getEnergy();
    std::cout<<"Energy in iteration "<<iter<<": "<<E<<"\n";
    std::cout<<"rho: "<<ro<<" r: "<<rr<<" s "<<ss<<"\n";

    //r = Ax+Bz-c primal residual
    //s = ro A^TB(z^{k+1} - z^k)  dual residual
    //a simple scheme that often workds well is
    //ro^{k+1} = tau * ro^k  if ||r^k||>mu||s^k||
    //ro^{k+1} = ro^k/tau  if   ||s^k||>mu||r^k||
    //typical choice mu=10,  tau = 2.
    //here r = \sum(x-z)^2 , s = \sum C_k(z^{k+1}-z^k)^2, where C_z is the
    //count of how many times z_k appears in constraints.
    float mu = 10, tau= 2;
//    ss *= ro;
//    if(rr > mu*ss/1000 || std::isnan(E)){
//      ro *= tau;
//      for(size_t ii = 0;ii<l.size();ii++){
//        Element * ele = e->elements[ii];
//        for(size_t jj = 0;jj<ele->getNumNodes();jj++){
//          l[ii][jj] /= tau;
//        }
//      }
//      for(auto ii = 0;ii<cl.size();ii++){
//        cl[ii] /= tau;
//      }
//    }
//    if(mu*rr < ss/1000){
//      ro /= tau;
//      for(size_t ii = 0;ii<l.size();ii++){
//        Element * ele = e->elements[ii];
//        for(size_t jj = 0;jj<ele->getNumNodes();jj++){
//          l[ii][jj] *= tau;
//        }
//      }
//    }

    if(std::abs(prevE-E) < tol){
      break;
    }
    prevE = E;

    tt = clock();
    out<<(tt-tt0)/(CLOCKS_PER_SEC/1000.0);
    float ene = e->getEnergy();
    out<<" "<<ene<<"\n";

  }
}

void ADMMStepper::Step(World * world)
{

  for(size_t ii = 0;ii<world->element.size();ii++){
    Step(world->element[ii]);

  }
}

ADMMStepper::ADMMStepper():ro(1),tol(0.001),xtol(0.001)
{
}
