#include "AdmmNoSpring.hpp"
#include "Element.hpp"
#include "ElementHex.hpp"
#include "ElementMesh.hpp"
#include "MaterialQuad.hpp"
#include "StrainLin.hpp"
#include "MatrixXd.hpp"
#include "LinSolve.hpp"
#include "ArrayUtil.hpp"
#include "femError.hpp"
#include "Quadrature.hpp"

#include <time.h>

float AdmmNoSpring::getEnergy(ElementMesh * eMesh, int eIdx)
{
  Element * e = eMesh->e[eIdx];
  float E = eMesh->getEnergy(eIdx);
  Eigen::VectorXf UZ = Eigen::VectorXf::Zero(3 * e->nV());
  for (int ii = 0; ii<e->nV(); ii++){
    Vector3f diff = eMesh->x[e->at(ii)] - Z[e->at(ii)];
    for (int jj = 0; jj < 3; jj++){
      UZ[3*ii+jj] = diff[jj];
    }
    E += Vector3f::dot(y[eIdx][ii], diff);
  }
  E += UZ.transpose() * T[eIdx] * UZ;
  return E;
}

std::vector<Vector3f>
AdmmNoSpring::getForces(ElementMesh * eMesh, int eIdx)
{
  ElementHex * ele = (ElementHex *)eMesh->e[eIdx];
  std::vector<Vector3f> ff = eMesh->getForce(eIdx);
  Eigen::VectorXf Uz = Eigen::VectorXf::Zero(3 * ele->nV());

  for (int ii = 0; ii < ele->nV(); ii++){
    for (int jj = 0; jj < 3; jj++){
      Uz[3 * ii + jj] = Z[ele->at(ii)][jj] - eMesh->x[ele->at(ii)][jj];
    }
  }
  //continuity force
  Eigen::VectorXf fCont = T[eIdx]*Uz;
  for (unsigned int ii = 0; ii<ff.size(); ii++){
    for (int jj = 0; jj < 3; jj++){
      ff[ii][jj] += fCont[3 * ii + jj];
    }
    //lagrange multiplier force
    ff[ii] -= y[eIdx][ii];
  }
  return ff;
}

MatrixXd
AdmmNoSpring::stiffness(ElementMesh *mesh, int eIdx)
{
  MatrixXd K = mesh->getStiffness(eIdx);
  
  //std::ofstream out("debug.txt");
  //K.print(out);
  //out << "\n";
  //out.close();

  for (int ii = 0; ii<K.mm; ii++){
    for (int jj = 0; jj < K.nn; jj++){
      K(ii, jj) += T[eIdx](ii, jj);
    }
  }
  return K;
}

void AdmmNoSpring::minimizeElement(ElementMesh * m, Element * ele,
  int eIdx)
{
  float E = 0;
  float h = 1;
  int NSteps = 1;
  int ndof = 3 * ele->nV();

  //std::ofstream out("debug.txt");
  for (int iter = 0; iter<NSteps; iter++){
    std::vector<Vector3f> force = getForces(m, eIdx);
    MatrixXd K = stiffness(m, eIdx);
    
    //K.print(out);
    //out << "\n";
    //out<<T<<"\n\n";

    for (int ii = 0; ii<ele->nV(); ii++){
      int vidx = ele->at(ii);
      if (m->fixed[vidx]){
        force[ii] = Vector3f(0, 0, 0);
      }
      for (int jj = 0; jj<3; jj++){
        int row = 3 * ii + jj;
        //K(row, row) += 100;
        bb[row] = force[ii][jj];
        if (m->fixed[vidx]){
          for (int kk = 0; kk<ndof; kk++){
            K(row, kk) = 0;
            K(row, row) = 10;
          }
        }
      }
    }

    //K.print(out);
    //out << "\n";
    //for (int ii = 0; ii < 3 * ele->nV(); ii++){
    //  out << bb[ii] << "\n";
    //}
    //out << "\n";
    
    linSolve(K, bb);
    
    //for (int ii = 0; ii < 3 * ele->nV(); ii++){
    //  out << bb[ii] << "\n";
    //}
    //out << "\n"; 
    //out.close();

    for (int ii = 0; ii<ele->nV(); ii++){
      for (int jj = 0; jj<3; jj++){
        force[ii][jj] = (float)bb[3 * ii + jj];
      }
    }

    float totalMag = 0;
    for (unsigned int ii = 0; ii<force.size(); ii++){
      int vidx = ele->at(ii);
      for (int jj = 0; jj<3; jj++){
        totalMag += std::abs(force[ii][jj]);
      }
    }
    if (totalMag<xtol){
      return;
    }

    float E = getEnergy(m, eIdx);

    //line search
    std::vector<Vector3f> x0, x;
    getEleX(eIdx, m, x0);
    float E1=E;
    while (1){
      x = x0;
      addmul(x, h, force);
      setEleX(eIdx, m, x);
      E1 = getEnergy(m, eIdx);
      // std::cout<<h<<" "<<E1<<"\n";
      if (E1>E || fem_error){
        fem_error = 0;
        h = 0.5f* h;
      }
      else{
        h = 1.1f*h;
        break;
      }
    }
  }
}

void AdmmNoSpring::init(ElementMesh *_m)
{
  m = _m;
  out.open("converge.txt");
  
  m->u = &u;

  float prevE = m->getEnergy();


  if (bb != 0){
    delete[] bb;
  }

  const Quadrature & q2d = Quadrature::Gauss2_2D;

  T.resize(m->e.size());
  for (unsigned int ee = 0; ee < m->e.size(); ee++){
    ElementHex * ele = (ElementHex*)m->e[ee];
    MaterialQuad * mat = (MaterialQuad*)m->m[m->me[ee]];
    StrainLin * ene = (StrainLin*)mat->e[0];
    Eigen::MatrixXf E = ene->EMatrix();   
    T[ee] = Eigen::MatrixXf::Zero(3 * ele->nV(), 3 * ele->nV());
    for (int ii = 0; ii < ele->nF(); ii++){
      Eigen::MatrixXf Tf = Eigen::MatrixXf::Zero(3 * ele->nV(), 3 * ele->nV());
      Eigen::MatrixXf N = ele->NMatrix(ii);
      for (unsigned int jj = 0; jj < q2d.x.size(); jj++){
        Vector3f quadp = ele->facePt(ii, q2d.x[jj]);
        Eigen::MatrixXf H = ele->HMatrix(quadp);
        //Eigen::MatrixXf B0 = ele->BMatrix(quadp, e->X);
        //Eigen::MatrixXf B = Eigen::MatrixXf::Zero(6, 3 * ele->nV());
        //only add contributions from the face
        //for (int kk = 0; kk < 4; kk++){
        //  int vidx = faces[ii][kk];
        //  B.block(0, 3 * vidx, 6, 3) = B0.block(0, 3 * vidx, 6, 3);
        //}
        //Tf += q2d.w[jj] * B.transpose() * E * B;
        Tf += q2d.w[jj] * H.transpose() * N * E * N.transpose() * H;
      }
      T[ee] += Tf;
    }
    float sideLen = m->X[ele->at(7)][0] - m->X[ele->at(0)][0];
    float area = sideLen*sideLen;
    T[ee] = area*T[ee];
  }

  bb = new double[3 * m->x.size()];
  u.resize(m->e.size());
  y.resize(u.size());
  for (size_t ii = 0; ii<m->e.size(); ii++){
    Element * ele = m->e[ii];
    u[ii].resize(ele->nV());
    y[ii].resize(u[ii].size());
    for (int jj = 0; jj<ele->nV(); jj++){
      u[ii][jj] = m->X[ele->at(jj)];
    }
  }
  Z = m->x;
}

int AdmmNoSpring::oneStep()
{
  clock_t tt, tt0;
  tt0 = clock();
  std::vector<Vector3f> x0 = m->x;

    //Z in the previous iteraiton.
    std::vector<Vector3f> Zk_1 = Z;

    //update u locally
    for (unsigned int ii = 0; ii<m->e.size(); ii++){
      Element * ele = m->e[ii];
      setEleX(ii, m, u[ii]);
      minimizeElement(m, ele, ii);
      getEleX(ii, m, u[ii]);
    }

    //restore x to consensus
    m->x = Z;

    std::vector<Eigen::Matrix3f> roSum(Z.size(), Eigen::Matrix3f::Zero());
    //update z closed form
    for (unsigned int ii = 0; ii<m->X.size(); ii++){
      Z[ii] = Vector3f(0, 0, 0);
    }

    //add per element variables and multipliers
    for (unsigned int ii = 0; ii<u.size(); ii++){
      Element * ele = m->e[ii];
      for (int jj = 0; jj<ele->nV(); jj++){
        int vIdx = ele->at(jj);
        Vector3f disp = u[ii][jj] - m->X[ele->at(jj)];
        Eigen::Vector3f uLocal(disp[0], disp[1], disp[2]);
        Eigen::Matrix3f Tblock = T[ii].block(3 * jj, 3 * jj, 3, 3);
        Eigen::Vector3f fEigen = Tblock * uLocal;
        Vector3f ff(fEigen[0], fEigen[1], fEigen[2]);
        Z[vIdx] +=  ff + y[ii][jj];
        roSum[vIdx] += Tblock;
      }
    }

    for (unsigned int ii = 0; ii< m->fe.size(); ii++){
      Vector3f force = m->fe[ii];
      Z[ii] += force;
    }
    
    for (unsigned int ii = 0; ii < Z.size(); ii++){
      Eigen::Vector3f ez(Z[ii][0], Z[ii][1], Z[ii][2]);
      ez = roSum[ii].colPivHouseholderQr().solve(ez);
      Z[ii] = Vector3f(ez[0], ez[1],ez[2]) + m->X[ii];
    }

    //fix constraints
    for (unsigned int ii = 0; ii<m->fixed.size(); ii++){
      if (m->fixed[ii]){
        Z[ii] = x0[ii];
      }
    }

    //line search for best delta Z magnitude
    float E = m->getEnergy();
    
    std::cout << "Energy: " << E << "\n";
    
    float ene1 = E;
    for (unsigned int ii = 0; ii<Z.size(); ii++){
      Z[ii] = Z[ii] - Zk_1[ii];
    }
    float hh = 1.0f;
    while (1){
      m->x = Zk_1;
      addmul(m->x, hh, Z);
      ene1 = m->getEnergy();
      if (ene1<E && fem_error == 0){
        Z = m->x;
        break;
      }
      else{
        hh = hh / 2;
        if (hh<1e-15){
          m->x = Zk_1;
          Z = Zk_1;
          break;
        }
      }
    }
    std::cout << "hh " << hh << "\n";
    if (prevE - ene1 < tol){
  //    break;
    }
    prevE = ene1;
    
    //update multiplier for elements
    for (unsigned int ii = 0; ii<u.size(); ii++){
      Element * ele = m->e[ii];
      Eigen::VectorXf UZ = Eigen::VectorXf::Zero(3 * ele->nV());
      for (int jj = 0; jj < ele->nV(); jj++){
        Vector3f diff = u[ii][jj] - Z[ele->at(jj)];
        for (int kk = 0; kk < 3; kk++){
          UZ[3 * jj + kk] = diff[kk];
        }
      }
      UZ = T[ii]*UZ;
      for (int jj = 0; jj<ele->nV(); jj++){
        for (int kk = 0; kk < 3; kk++){
          y[ii][jj][kk] += UZ[3 * jj + kk];
        }
      }
    }

    tt = clock();
    out << (tt - tt0) / (CLOCKS_PER_SEC / 1000.0);
    return 0;
}

AdmmNoSpring::~AdmmNoSpring()
{
  if (bb != 0){
    delete bb;
  }
}

AdmmNoSpring::AdmmNoSpring() :bb(0)
, tol(0.001f), xtol(0.001f)
{
}
