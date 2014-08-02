/*
 * Copyright 1993-2014 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

/*
 * This sample implements a conjugate graident solver on GPU
 * using CUBLAS and CUSPARSE
 *
 */

// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Using updated (v2) interfaces to cublas and cusparse */
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include "helper_cuda.h"
#include "ConjugateGradientCuda.hpp"

ConjugateGradientCuda::ConjugateGradientCuda():tol(1e-5f),max_iter(100),devID(0),
  cublasHandle(0), cusparseHandle(0),descr(0),N(0),nz(0),
  d_col(0), d_row(0), d_val(0), d_x(0), d_r(0), d_p(0), d_Ax(0)
{}

int ConjugateGradientCuda::initCuda(int _N, int _nz, int * _I, int * _J)
{
  N=_N;
  nz=_nz;
  I=_I;
  J=_J;
  // This will pick the best possible CUDA capable device
  devID = findCudaDevice();
  //cudaDeviceProp deviceProp;
  //    checkCudaErrors(cudaGetDeviceProperties(&deviceProp, devID));
  if (devID < 0)
  {
    printf("no cuda gpu\n");
    return -1;
  }

  /* Get handle to the CUBLAS context */
  cublasStatus_t cublasStatus;
  cublasStatus = cublasCreate(&cublasHandle);

  checkCudaErrors(cublasStatus);

  /* Get handle to the CUSPARSE context */
  
  cusparseStatus_t cusparseStatus;
  cusparseStatus = cusparseCreate(&cusparseHandle);

  checkCudaErrors(cusparseStatus);

  cusparseStatus = cusparseCreateMatDescr(&descr);

  checkCudaErrors(cusparseStatus);

  cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);

  checkCudaErrors(cudaMalloc((void **)&d_col, nz*sizeof(int)));
  checkCudaErrors(cudaMalloc((void **)&d_row, (N+1)*sizeof(int)));
  checkCudaErrors(cudaMalloc((void **)&d_val, nz*sizeof(float)));
  checkCudaErrors(cudaMalloc((void **)&d_x, N*sizeof(float)));
  checkCudaErrors(cudaMalloc((void **)&d_r, N*sizeof(float)));
  checkCudaErrors(cudaMalloc((void **)&d_p, N*sizeof(float)));
  checkCudaErrors(cudaMalloc((void **)&d_Ax, N*sizeof(float)));

  cudaMemcpy(d_col, J, nz*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_row, I, (N+1)*sizeof(int), cudaMemcpyHostToDevice);
  return 0;
}

int ConjugateGradientCuda::solve(float * val, float * x, float * rhs)
{
    float a, b, na, r0, r1,dot;
    int k;
    float alpha, beta, alpham1;

    cublasStatus_t cublasStatus;
    
    cudaMemcpy(d_val, val, nz*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x, N*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, rhs, N*sizeof(float), cudaMemcpyHostToDevice);

    alpha = 1.0;
    alpham1 = -1.0;
    beta = 0.0;
    r0 = 0.;

    cusparseScsrmv(cusparseHandle,CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_x, &beta, d_Ax);

    cublasSaxpy(cublasHandle, N, &alpham1, d_Ax, 1, d_r, 1);
    cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);

    k = 1;

    while (r1 > tol*tol && k <= max_iter)
    {
        if (k > 1)
        {
            b = r1 / r0;
            cublasStatus = cublasSscal(cublasHandle, N, &b, d_p, 1);
            cublasStatus = cublasSaxpy(cublasHandle, N, &alpha, d_r, 1, d_p, 1);
        }
        else
        {
            cublasStatus = cublasScopy(cublasHandle, N, d_r, 1, d_p, 1);
        }

        cusparseScsrmv(cusparseHandle, CUSPARSE_OPERATION_NON_TRANSPOSE, N, N, nz, &alpha, descr, d_val, d_row, d_col, d_p, &beta, d_Ax);
        cublasStatus = cublasSdot(cublasHandle, N, d_p, 1, d_Ax, 1, &dot);
        a = r1 / dot;

        cublasStatus = cublasSaxpy(cublasHandle, N, &a, d_p, 1, d_x, 1);
        na = -a;
        cublasStatus = cublasSaxpy(cublasHandle, N, &na, d_Ax, 1, d_r, 1);

        r0 = r1;
        cublasStatus = cublasSdot(cublasHandle, N, d_r, 1, d_r, 1, &r1);
        cudaThreadSynchronize();
        printf("iteration = %3d, residual = %e\n", k, sqrt(r1));
        k++;
    }

    cudaMemcpy(x, d_x, N*sizeof(float), cudaMemcpyDeviceToHost);

    float rsum, diff, err = 0.0;

    for (int i = 0; i < N; i++)
    {
        rsum = 0.0;

        for (int j = I[i]; j < I[i+1]; j++)
        {
            rsum += val[j]*x[J[j]];
        }

        diff = fabs(rsum - rhs[i]);

        if (diff > err)
        {
            err = diff;
        }
    }

//    printf("Test Summary:  Error amount = %f\n", err);
    return 0;
}

int ConjugateGradientCuda::clearCuda()
{
  cusparseDestroy(cusparseHandle);
  cublasDestroy(cublasHandle);

  cudaFree(d_col);
  cudaFree(d_row);
  cudaFree(d_val);
  cudaFree(d_x);
  cudaFree(d_r);
  cudaFree(d_p);
  cudaFree(d_Ax);

  // cudaDeviceReset causes the driver to clean up all state. While
  // not mandatory in normal operation, it is good practice.  It is also
  // needed to ensure correct operation when the application is being
  // profiled. Calling cudaDeviceReset causes all profile data to be
  // flushed before the application exits
  cudaDeviceReset();
  return 0;
}
