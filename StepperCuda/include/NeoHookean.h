#ifndef NEOHOOKEAN_H
#define NEOHOOKEAN_H
#include "cuda_runtime.h"
#include "float3x3.h"
#include "float3x3_helper.h"
#include "math_functions.h"
#include "cuPrintf.cu"
__device__ float param[2] = {34482,310344};

__device__ float FTFTrace(const float3x3 & F)
{
	float r = 0;
	for(int ii = 0;ii<9;ii++){
		r += F.m[ii] * F.m[ii];
	}
	return r;
}

///@return -1 if inverted. Energy is assumed to be non-negative.
__device__ float StrainEnergy(const float3x3 & F)
{
	float det = F.det();
	if(det<1e-10){
		return -1;
	}
	float I1 = FTFTrace(F);
	float JJ = __logf(det);
	float mu = param[0],lambda=param[1];
	float Psi = (mu/2) * (I1-3) - mu*JJ + (lambda/2)*JJ*JJ;
	return Psi;
}

__device__ void PK1(const float3x3 & F , float3x3 & PP)
{
	float JJ = __logf(F.det());
	float3x3 Finv;
	F.inverse(Finv);
	Finv.transpose();
	float mu = param[0],lambda=param[1];
	PP = mu*(F-Finv) + lambda*JJ*Finv;
}

__device__ void dPdx(const float3x3 & F, const float3x3 & dF, float3x3 & dP)
{
	dP.setZero();
	float JJ = __logf(F.det());
	float3x3 Finv;
	F.inverse(Finv);
	float trace = (Finv*dF).trace();
	Finv.transpose();
	float3x3 dFT = dF;
	dFT.transpose();
	float mu = param[0],lambda=param[1];
	dP = mu*dF;
	float c1 = mu-lambda * JJ;
	dP += c1 * Finv * dFT * Finv;
	dP += lambda * trace * Finv;
}

#endif