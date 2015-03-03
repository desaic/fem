#ifndef HEXELE_H
#define HEXELE_H
#include "cuda_runtime.h"
#include "vector_types.h"

#define MAX_QUADRATURE 64
#define NVERT 8
#define KSIZE (3*NVERT)

typedef struct {
	float3 P[MAX_QUADRATURE];
	float W[MAX_QUADRATURE];
	int N;
}QuadratureCuda;

struct KMat
{
	float K[KSIZE][KSIZE];
};

struct ADMMInfo
{
	float ro;
	float3 Z[NVERT];
  float3 y[NVERT];
};

__global__ void HexTest(float3 * xx);
__global__ void HexStiffnessTest(const float3 * XX, const float3 * xx, KMat * Ks,
								 float3 * ff, float * E);
__global__ void admmMinTest(const float3 * XX, float3 * xx, const ADMMInfo * admm);
///@brief each element keeps a duplicated copy of Z coord stored in admm
__global__ void admmMinEleDup(const float3 * XX, float3 * xx, const ADMMInfo * admm);

__global__ void batchCG(const float * mats, const float * bs, float * xs, int nmat);

__global__ void GetEnergy(const float3 * XX,  const ADMMInfo * admm, float * E);
///@brief compute internal forces
__global__ void GetIntForce(const float3 * XX, const ADMMInfo * admm, float3 * ff);

void initHexEle();
#endif
