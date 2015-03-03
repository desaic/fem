#ifndef CGCUDA_H
#define CGCUDA_H
#include "cuda_runtime.h"

#define MATSIZE 24

__device__ void addMult(const float * a, const float * b, float c, float * r)
{
	for(int ii = 0; ii<MATSIZE; ii++){
		r[ii] = a[ii] + c * b[ii];
	}
}

__device__ float dot(const float * a, const float * b)
{
	float dd=0;
	for(int ii = 0; ii<MATSIZE; ii++){
		dd+=a[ii] * b[ii];
	}
	return dd;
}

__device__ void mulMatVec(const float * A,const float * x, float * r)
{
	for (int ii = 0; ii < MATSIZE; ii++)
	{
		float prodEntry = 0.0f;
		for (int jj = 0; jj < MATSIZE; jj++)
		{
			prodEntry += A[jj + ii*MATSIZE] * x[jj];
		}
		r[ii] = prodEntry;
	}
}

///@param N square matrix dimension.
__device__ void CG(const float * mat, const float * bb, float * xx)
{
	///@TODO need better solver
	const int maxIter=15;
	const float tol=0.001;
	float rsold=0;
	//conjugate gradient dir and residual vector
	float pp[MATSIZE], rr[MATSIZE], Ap[MATSIZE];
	//initial guess 0.
	for(int ii =0;ii<MATSIZE;ii++){
		xx[ii] = 0;
	}
	for(int ii =0;ii<MATSIZE;ii++){
		pp[ii] = bb[ii];
		rsold += pp[ii] * pp[ii];
	}
	for(int ii =0;ii<MATSIZE;ii++){
		rr[ii] = pp[ii];
	}
	if(rsold<tol){
		return;
	}
	for (int ii=0;ii<maxIter;ii++){
		//rsnew = A*p temporarily
		mulMatVec(mat, pp, Ap);
		float alpha = rsold/dot(pp, Ap);
		addMult(xx, pp,  alpha, xx);
		addMult(rr, Ap, -alpha, rr);
		float rsnew = dot(rr, rr);
		if(rsnew<tol){
			break;
		}
		addMult(rr, pp, rsnew/rsold,pp);
		rsold = rsnew;
	}
}

__global__ void batchCG(const float * mats, const float * bs, float * xs, int nmat)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	int NN = MATSIZE;
	float mat[MATSIZE*MATSIZE];
	float bb[MATSIZE];
	float xx[MATSIZE];
	//assemble sample SPD matrix.
	for(int jj = 0;jj<NN;jj++){
		for(int kk = 0;kk<NN;kk++){
			mat[jj*NN+kk]=0.01;
		}
		int kk=jj;
		mat[jj*NN+kk] = 2.1f;
		if(jj>0){
			mat[jj*NN+kk-1]=-1;
		}
		if(jj<NN-1){
			mat[jj*NN+kk+1]=-1;
		}
		bb[jj] = 1;
		xx[jj] = 0;
	}
	CG(mat, bb, xx);
	for(int ii = 0;ii<NN;ii++){
		xs[i*NN + ii] = xx[ii];
	}
}
#endif