#include "UnitTests.hpp"
#include "Element.hpp"
#include "ElementRegGrid.hpp"
#include "ElementHex.hpp"
#include "MaterialQuad.hpp"
#include "StrainEneNeo.hpp"
#include "StrainLin.hpp"
#include "StrainCorotLin.hpp"
#include "Quadrature.hpp"
#include "StepperNewton.hpp"
//#include "LinSolveCusp.hpp"
#include "Timer.hpp"

#include <Eigen/Dense>
#include <iostream>
#include <fstream>

using namespace Eigen;

void forceTestHelper(StrainEne * ene);
void stiffnessTestHelper(StrainEne * ene);

void cudaLinTest()
{
  int nx = 1, ny = 1, nz = 1;
  ElementRegGrid * grid = new ElementRegGrid(nx, ny, nz);
  StrainEne * ene = 0;
  ene = new StrainLin();
  ene->param[0] = 10000;
  ene->param[1] = 100000;
  MaterialQuad * material = new MaterialQuad(ene);
  grid->m.push_back(material);
  
  std::ofstream out("debug.txt");

  Vector3f ff(10, -50, 0);
  ff = (1.0f / (nx*nz)) * ff;

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
  int status = grid->check();
  std::cout << "Check Mesh: " << status << "\n";
  std::vector<int> I, J;
  std::vector<float> val;
  grid->stiffnessPattern(I, J);
  grid->getStiffnessSparse(val);
  float * x = new float[I.size()-1];
  std::vector<Vector3f> f = grid->getForce();
  for (unsigned int row = 0; row < I.size() - 1; row++){
    int nCol = I[row + 1] - I[row];
    int vi = row / 3;
    for (int jj = 0; jj < nCol; jj++){
      int idx = I[row] + jj;
      int col = J[idx];
      int vj = col / 3;
      
      if (grid->fixed[vi] || grid->fixed[vj]){
        val[idx] = 0;
        if (row == col){
          val[idx] = 100;
        }
      }
      out << val[idx] << " ";
    }
    out << "\n";
  }
  for (unsigned int ii = 0; ii < f.size(); ii++){
    if (grid->fixed[ii]){
      f[ii] = Vector3f::Zero();
    }
    for (int jj = 0; jj < 3; jj++){
      x[3 * ii + jj] = f[ii][jj];
      out << x[3 * ii + jj] << "\n";
    }
  }
  out.close();
//  LinSolveCusp solver;
//  solver.init(I, J);
//  Timer timer;
//  timer.start();
//  solver.solve(val, x);
//  for (unsigned int ii = 0; ii < f.size(); ii++){
//    for (int jj = 0; jj < 3; jj++){
//      f[ii][jj] = x[3 * ii + jj];
//    }
//    std::cout << f[ii][0] << " " << f[ii][1] << " " << f[ii][2] << "\n";
//  }
//  timer.end();
//  std::cout << timer.getSeconds();
  delete[]x;
}

void stiffnessTest(int matModel)
{
  StrainEne * ene = 0;
  switch (matModel){
  case 0:
    ene = new StrainLin();
    ene->param[0] = 34482;
    ene->param[1] = 310344;
    break;
  case 1:
    ene = new StrainCorotLin();
    ene->param[0] = 34482;
    ene->param[1] = 310344;
    break;
  case 2:
    ene = new StrainEneNeo();
    ene->param[0] = 34482;
    ene->param[1] = 310344;
    break;
  default:
    return;
  }
  stiffnessTestHelper(ene);
  delete ene;

}

void stiffnessTestHelper(StrainEne * ene)
{
  float h = 0.001f;
  ElementRegGrid * grid = new ElementRegGrid(1, 1, 1);
  MaterialQuad * material = new MaterialQuad(ene);
  material->init(grid);
  grid->m.push_back(material);
 // grid->x[1][0] += 0.1f;
  //grid->x[3][1] += 0.2f;
  MatrixXf KAna = grid->getStiffness();
  int nVar = (int)grid->X.size();
  MatrixXf K(3*nVar,3*nVar);
  //check each partial derivative
  for(unsigned int ii = 0;ii<grid->x.size();ii++){
    for(int jj = 0; jj<3; jj++){
      grid->x[ii][jj] -= h;
      std::vector<Vector3f>fminus = grid->getForce();
      grid->x[ii][jj] += 2*h;
      std::vector<Vector3f>fplus = grid->getForce();
      grid->x[ii][jj] -=h;
      for(unsigned int kk = 0;kk<fminus.size();kk++){
        fplus[kk] -= fminus[kk];
        fplus[kk] /= 2*h;
        for(int ll = 0;ll<3;ll++){
          K(3*kk+ll,3*ii+jj)=-fplus[kk][ll];
        }
      }
    }
  }

  std::cout<<"Ana K:\n";
  std::cout << KAna << "\nNum K:\n" << K << "\n";
  float maxErr = 0;
  for(int ii = 0;ii<K.rows();ii++){
    for(int jj =0 ;jj<K.cols();jj++){
      float err = (float)std::abs(KAna(ii,jj)-K(ii,jj));
      if(err>maxErr){
        maxErr = err;
      }
    }
  }
  std::cout<<"max err "<<maxErr<<"\n";

  delete grid;
  delete material;
}

void linearTMatrixTest(StrainLin * ene)
{
  ElementRegGrid * grid = new ElementRegGrid(1, 1, 1);
  MaterialQuad * material = new MaterialQuad(ene);
  grid->m.push_back(material);
  grid->x[1][0] += 0.1f;
  grid->x[3][1] += 0.2f;
  MatrixXf K = grid->getStiffness();

  //linear material stiffness
  ElementHex * ele = (ElementHex*)grid->e[0];
  const Quadrature & q = Quadrature::Gauss2;
  Eigen::MatrixXf Ka = Eigen::MatrixXf::Zero(3 * ele->nV(), 3 * ele->nV());
  Eigen::MatrixXf E = ene->EMatrix();
  Eigen::VectorXf U = Eigen::VectorXf::Zero(3 * ele->nV());

  for (int ii = 0; ii < ele->nV(); ii++){
    for (int jj = 0; jj < 3; jj++){
      U[3 * ii + jj] = grid->x[ii][jj] - grid->X[ii][jj];
    }
  }

  for (unsigned int ii = 0; ii < q.x.size(); ii++){
    Eigen::MatrixXf B = ele->BMatrix(q.x[ii], grid->X);
    Eigen::MatrixXf ss = E*B*U;
    //std::cout <<"sigma:\n"<< ss << "\n";

    Matrix3f F = ele->defGrad(q.x[ii], grid->X, grid->x);
    Matrix3f P = ene->getPK1(F);
    //std::cout << "P:\n";
    //P.print();

    Ka += q.w[ii] * B.transpose() * E * B;
    //std::cout << "B:\n" << B << "\n";
  }

  //std::cout << "E:\n" << E << "\n";
  //std::cout << "alt K:\n";
  //std::cout << Ka << "\n";
  float maxErr = 0;
  for (int ii = 0; ii<K.rows(); ii++){
    for (int jj = 0; jj<K.cols(); jj++){
      float err = (float)std::abs(Ka(ii, jj) - K(ii, jj));
      if (err>maxErr){
        maxErr = err;
      }
    }
  }

  std::cout << "max err " << maxErr << "\n";

  //assemble boundary matrix HNEB
  std::ofstream out("T.txt");

  // 2 point quadrature is accurate enough
  const Quadrature & q2d = Quadrature::Gauss2_2D;
  Eigen::MatrixXf T = Eigen::MatrixXf::Zero(3 * ele->nV(), 3 * ele->nV());
  for (int ii = 0; ii < ele->nF(); ii++){
    Eigen::MatrixXf Tf = Eigen::MatrixXf::Zero(3 * ele->nV(), 3 * ele->nV());
    Eigen::MatrixXf N = ele->NMatrix(ii);
    N.block(0, 3, 3, 3) = Eigen::MatrixXf::Zero(3, 3);
    //std::cout << "N:\n"<<N << "\n";
    for (unsigned int jj = 0; jj < q2d.x.size(); jj++){
      Vector3f quadp = ele->facePt(ii, q2d.x[jj]);
      Eigen::MatrixXf B0 = ele->BMatrix(quadp, grid->X);
      Eigen::MatrixXf B = Eigen::MatrixXf::Zero(6, 3 * ele->nV());
      //only add contributions from the face
      for (int kk = 0; kk < 4; kk++){
        int vidx = faces[ii][kk];
        B.block(0, 3 * vidx, 6, 3) = B0.block(0, 3 * vidx, 6, 3);
      }
      //B=B0;
      Eigen::MatrixXf H = ele->HMatrix(quadp);
      //std::cout << "H:\n" << H.transpose() << "\n";
      //std::cout << "B:\n" << B.transpose() << "\n";
      //std::cout << "traction:\n";
      //Tf += q2d.w[jj] * H.transpose() * N * E * B;
      Tf += q2d.w[jj] * H.transpose() * N * E * N.transpose() * H;
      //Tf += q2d.w[jj] * B.transpose() * E * B;
      Eigen::Vector3f surfF = (N * E * B * U);
      //std::cout << surfF << "\n";
      Matrix3f F = ele->defGrad(quadp, grid->X, grid->x);
      Matrix3f P = ene->getPK1(F);
      Vector3f surfF1 = P * Vector3f(facen[ii][0], facen[ii][1], facen[ii][2]);
      std::cout << surfF1[0] << " " << surfF1[1] << " " << surfF1[2] << "\n";
    }
    //out << Tf << "\n\n";
    T += Tf;
  }
  //out << T << "\n\n";
  //out << Ka << "\n";
  out.close();
}

void forceTest(int matModel)
{
  StrainEne * ene = 0;
  switch (matModel){
  case 0:
    ene = new StrainLin();
    ene->param[0] = 34482;
    ene->param[1] = 310344;
    break;
  case 1:
    ene = new StrainCorotLin();
    ene->param[0] = 34482;
    ene->param[1] = 310344;
    break;
  case 2:
    ene = new StrainEneNeo();
    ene->param[0] = 34482;
    ene->param[1] = 310344;
    break;
  default:
    return;
  }
  forceTestHelper(ene);
  delete ene;
}

void forceTestHelper(StrainEne * ene)
{
  int nx = 1,ny=1,nz=1;
  ElementRegGrid * em = new ElementRegGrid(nx,ny,nz);
  MaterialQuad material(ene);
  em->m.push_back(&material);
  
  for(int ii = 0;ii<nx;ii++){
    for(int jj =0;jj<nz;jj++){
      int eidx= em->GetEleInd(ii,0,jj);
      int aa[4] = {0,1,4,5};
      int bb[4] = {2,3,6,7};
      for(int kk = 0;kk<4;kk++){
        int vidx =em->e[eidx]->at(aa[kk]);
        //      em->fixed[vidx] = 1;
        vidx = em->e[eidx]->at(bb[kk]);
        //    em->fe[vidx] = ff;
      }
    }
  }
  em->check();
  material.init(em);
  em->x[0][0] += 0.2f;
  em->x[1][2] -= 0.3f;
  float h = 0.0001f;

  std::vector<Vector3f>force = em->getForce();
  for(size_t ii = 0;ii<em->x.size();ii++){
    for(int jj = 0; jj<3; jj++){
      em->x[ii][jj] -= h;
      double Eminus = em->getEnergy();
      em->x[ii][jj] += 2*h;
      double Eplus = em->getEnergy();
      std::cout<<"Energy diff: "<<Eplus-Eminus<<"\n";
      double numF = (Eplus - Eminus )/(2*h);
      std::cout<<"Analytic derivative:\n"<<-force[ii][jj];
      std::cout<<"\nCentral diff:\n"<<numF<<"\n--------\n";
    }
  }
  delete em;
}
