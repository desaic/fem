#include "device_launch_parameters.h"
#include "HexEle.h"
//#include "cuprintf.cuh"
#include "cuprintf.cu"
#include "float3x3.h"
#include "float3x3_helper.h"
#include "helper_math.h"
#include "NeoHookean.h"
#include "CGCuda.h"
#include "Quadrature.hpp"

__constant__ QuadratureCuda quadrature;
__constant__ float admmtol = 0.01;
///@brief error tolerance in energy for admm min ele.
__constant__ float etol = 1e-3;
__constant__ int nLinStep = 20;
__constant__ int sw[8][3] =
{ { -1, -1, -1 },
{ -1, -1, 1 },
{ -1, 1, -1 },
{ -1, 1, 1 },
{ 1, -1, -1 },
{ 1, -1, 1 },
{ 1, 1, -1 },
{ 1, 1, 1 }
};

__device__ void Stiffness(const float3 & p, const float3 *XX, const float3 * xx,
  float weight,
  float K[][KSIZE]);

///@param XX an array of 8 float3 s. Ordered z-> y ->x
__device__ float3 ShapeFunGrad(int ii, const float3 & xx,
  const float3 * XX)
{
  float3 size = 4 * (XX[7] - XX[0]);
  float3 grad;
  size.x = 1.0f / (size.x);
  size.y = 1.0f / (size.y);
  size.z = 1.0f / (size.z);

  grad.x = size.x * sw[ii][0] * (1 + sw[ii][1] * xx.y) * (1 + sw[ii][2] * xx.z);
  grad.y = size.y * sw[ii][1] * (1 + sw[ii][0] * xx.x) * (1 + sw[ii][2] * xx.z);
  grad.z = size.z * sw[ii][2] * (1 + sw[ii][0] * xx.x) * (1 + sw[ii][1] * xx.y);

  return grad;
}

__device__ void GetDefGrad(const float3 * XX, const float3 * xx, float3 p,
  float3x3& F)
{
  F.setIdentity();
  for (int ii = 0; ii<NVERT; ii++){
    float3 grad = ShapeFunGrad(ii, p, XX);
    float3x3 outer;
    outer.outerProd(xx[ii] - XX[ii], grad);
    F += outer;
  }
}

///@return -1 if inverted. Energy is assumed to be non-negative.
__device__ float GetEnergy(const float3 *XX, const float3 * xx)
{
  float energy = 0;
  for (int ii = 0; ii<quadrature.N; ii++){
    float3x3 F;
    GetDefGrad(XX, xx, quadrature.P[ii], F);
    float ee = StrainEnergy(F);
    if (ee <= -1){
      //for(int ii = 0;ii<9;ii++){
      //	cuPrintf("%.7f ", F.m[ii]);
      //}
      return -1;
    }
    energy += quadrature.W[ii] * ee;
  }
  float3 size = XX[7] - XX[0];
  float vol = size.x * size.y * size.z;
  return  vol*energy;
}

__device__ void ShapeFun(const float3 & p, float * ww)
{
  for (int ii = 0; ii<NVERT; ii++){
    ww[ii] = (1.0 / 8) * (1 + sw[ii][0] * p.x)
      *(1 + sw[ii][1] * p.y)
      *(1 + sw[ii][2] * p.z);
  }
}

__device__ void GetForce(const float3 *XX, const float3 * xx, float3 * ff)
{

  float3x3 F;
  float3 size = XX[7] - XX[0];
  float Nshape[NVERT];
  float vol = size.x * size.y * size.z;
  for (int ii = 0; ii<NVERT; ii++){
    ff[ii].x = 0;
    ff[ii].y = 0;
    ff[ii].z = 0;
  }
  for (int jj = 0; jj<quadrature.N; jj++){
    ShapeFun(quadrature.P[jj], Nshape);
    GetDefGrad(XX, xx, quadrature.P[jj], F);
    float3x3 PP;
    PK1(F, PP);
    for (int ii = 0; ii<NVERT; ii++){
      float3 gradN = ShapeFunGrad(ii, quadrature.P[jj], XX);
      ff[ii] -= vol* quadrature.W[jj] * PP * gradN;;
    }
  }
}

__device__ void Stiffness(const float3 * XX, const float3 * xx,
  float K[][KSIZE])
{
  float3 size = XX[7] - XX[0];
  float vol = size.x*size.y*size.z;

  for (int ii = 0; ii<KSIZE; ii++){
    for (int jj = 0; jj<KSIZE; jj++){
      K[ii][jj] = 0;
    }
  }
  for (int ii = 0; ii<quadrature.N; ii++){
    Stiffness(quadrature.P[ii], XX, xx, quadrature.W[ii] * vol, K);
  }
}

///@param p quadrature point in natural coordinates.
///@param weight Quadrature weight at point p.
///@TODO probably there is a more efficient calculation
__device__ void Stiffness(const float3 & p, const float3 *XX, const float3 * xx,
  float weight,
  float K[][KSIZE])
{
  float3 dN[8];
  float3x3 F;
  float3x3 dF;
  GetDefGrad(XX, xx, p, F);
  for (int ii = 0; ii<NVERT; ii++){
    dN[ii] = ShapeFunGrad(ii, p, XX);
  }
  for (int ii = 0; ii<8; ii++){
    for (int jj = 0; jj<3; jj++){
      //int col = 3*at[ii]+jj;
      int col = 3 * ii + jj;
      dF.setZero();
      dF.m[3 * jj] = dN[ii].x;
      dF.m[3 * jj + 1] = dN[ii].y;
      dF.m[3 * jj + 2] = dN[ii].z;
      float3x3 dP;
      dPdx(F, dF, dP);
      for (int vv = 0; vv<8; vv++){
        float3 dfdxi = dP*dN[vv];
        //K[3*at[vv]][col]   += weight*dfdxi.x;
        K[3 * vv][col] += weight*dfdxi.x;
        K[3 * vv + 1][col] += weight*dfdxi.y;
        K[3 * vv + 2][col] += weight*dfdxi.z;
      }
    }
  }
}

__device__ float admmEnergy(const float3 * XX, const float3 * xx,
  const ADMMInfo & admm)
{
  float E = GetEnergy(XX, xx);
  for (int ii = 0; ii<NVERT; ii++){
    float3 diff = xx[ii] - admm.Z[ii];
    E += 0.5 * admm.ro * dot(diff, diff);
    E += dot(admm.y[ii], diff);
  }
  return E;
}

__device__ void admmForce(const float3 * XX, const float3 * xx,
  const ADMMInfo & admm,
  float3 * ff)
{
  GetForce(XX, xx, ff);
  for (int ii = 0; ii<NVERT; ii++){
    ff[ii] -= admm.ro * (xx[ii] - admm.Z[ii]);
    ff[ii] -= admm.y[ii];
  }
}

__device__ void admmStiffness(const float3 * XX, const float3 * xx,
  const ADMMInfo & admm,
  float K[][KSIZE])
{
  Stiffness(XX, xx, K);
  for (int ii = 0; ii<KSIZE; ii++){
    K[ii][ii] += admm.ro;
  }
}

__device__ void arrayCpy(float3 * dst, const float3 * src, int N)
{
  for (int ii = 0; ii<N; ii++){
    dst[ii] = src[ii];
  }
}

__device__ void addMult(float3 * dst, const float3 * dx, float hh, int N)
{
  for (int ii = 0; ii<N; ii++){
    dst[ii] += hh*dx[ii];
  }
}

__device__ void admmMinEle(const float3 * XX, float3 * xx,
  const ADMMInfo & admm)
{
  float E = 0;
  float hh = 1;
  int NSteps = 2;
  float change = 0;
  float3 ff[NVERT];
  float3 dx[NVERT];
  float K[KSIZE][KSIZE];
  E = admmEnergy(XX, xx, admm);
  for (int ii = 0; ii<NSteps; ii++){
    //		cuPrintf("%d %.3f\n",ii, E);
    admmForce(XX, xx, admm, ff);
    for (int jj = 0; jj<NVERT; jj++){
      dx[jj].x = 0;
      dx[jj].y = 0;
      dx[jj].z = 0;
    }
    admmStiffness(XX, xx, admm, K);
    CG((const float*)K, (const float*)ff, (float*)dx);
    for (int jj = 0; jj<NVERT; jj++){
      change += abs(dx[jj].x);
      change += abs(dx[jj].y);
      change += abs(dx[jj].z);
    }
    //cuPrintf("f: %.4f\n", ff[0].x);
    //cuPrintf("dx %.4f\n", K[0][0]);
    //cuPrintf("dx %.4f\n", dx[0].x);
    //simple line search
    //ff=x0 x in previous iteration
    arrayCpy(ff, xx, NVERT);
    for (int jj = 0; jj<nLinStep; jj++){
      if (change * hh < etol){
        return;
      }
      addMult(xx, dx, hh, NVERT);
      float ene = admmEnergy(XX, xx, admm);
      //cuPrintf("%d %d\n", ii, jj);
      if (ene <= -1 || ene >= E){
        //cuPrintf("%.6f %0.6f\n", ene, E);
        //cuPrintf("%.4f\n", hh);
        hh *= 0.5;
        arrayCpy(xx, ff, NVERT);
      }
      else{
        if (abs(ene - E)<etol){
          return;
        }
        E = ene;
        break;
      }
    }
  }
}

__global__ void admmMinTest(const float3 * XX, float3 * xx, const ADMMInfo * admm)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  cuPrintfRestrict(0, 0);
  admmMinEle(XX, &(xx[NVERT*idx]), *admm);
}

__global__ void admmMinEleDup(const float3 * XX, float3 * xx, const ADMMInfo * admm)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  //cuPrintfRestrict(0,0);
  ADMMInfo admmLoc = admm[idx];
  admmMinEle(XX, &(xx[NVERT*idx]), admmLoc);
}

__global__ void GetEnergy(const float3 * XX, const ADMMInfo * admm, float * E)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  E[idx] = GetEnergy(XX, admm[idx].Z);
}

__global__ void GetIntForce(const float3 * XX, const ADMMInfo * admm, float3 * ff)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  GetForce(XX, admm[idx].Z, &(ff[NVERT*idx]));
}

__global__ void HexStiffnessTest(const float3 * XX, const float3 * xx, KMat * Ks,
  float3 * ff, float * E)
{
  cuPrintfRestrict(0, 0);
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  Stiffness(XX, xx, Ks[idx].K);
  E[idx] = GetEnergy(XX, xx);
  float3x3 F;
  float3x3 P;
  float3 p = make_float3(-0.5, -0.5, -0.5);
  GetDefGrad(XX, xx, p, F);
  GetForce(XX, xx, &(ff[8 * idx]));
  float3x3 dF;
  dF.setZero();
  float3 grad = ShapeFunGrad(6, p, XX);
  dF.m[6] = grad.x;
  dF.m[7] = grad.y;
  dF.m[8] = grad.z;
  dPdx(F, dF, P);
  for (int ii = 0; ii<3; ii++){
    ff[8 * idx + ii].x = P.m[3 * ii];
    ff[8 * idx + ii].y = P.m[3 * ii + 1];
    ff[8 * idx + ii].z = P.m[3 * ii + 2];

    //ff[8*idx + ii].x = F.m[3*ii];
    //ff[8*idx + ii].y = F.m[3*ii+1];
    //ff[8*idx + ii].z = F.m[3*ii+2];
  }
}

__global__ void HexTest(float3 * xx)
{
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  for (int ii = 0; ii<quadrature.N; ii++){
    xx[idx * 8 + ii].x = sw[ii][0];
  }
}

void initHexEle()
{
  QuadratureCuda qq;
  Quadrature quad=Quadrature::Gauss2;
  qq.N = 8;
  for (int ii = 0; ii<qq.N; ii++){
    qq.P[ii] = make_float3(quad.x[ii][0], quad.x[ii][1], quad.x[ii][2]);
    qq.W[ii] = quad.w[ii];
  }
  cudaMemcpyToSymbol(quadrature, &qq, sizeof(qq));
}
