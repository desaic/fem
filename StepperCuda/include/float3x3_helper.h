#ifndef FLOAT3X3_HELPER
#define FLOAT3X3_HELPER
#include "float3x3.h"
__device__ float3x3 operator*(float s, const float3x3 & mat)
{
	float3x3 b = mat;
	for(int ii = 0; ii<float3x3::MATLEN; ii++){
		b.m[ii] *= s;
	}
	return b;
}

__device__ float3x3 operator*(const float3x3 & a, const float3x3 & b)
{
	float3x3 c;
	c.setZero();
	for(int ii = 0;ii<3;ii++){
		for(int jj = 0;jj<3;jj++){
			for(int kk = 0;kk<3;kk++){
				c.m[3*ii+jj] += a.m[3*ii + kk] * b.m[3*kk + jj];
			}
		}
	}
	return c;
}

__device__ float3 operator*(const float3x3 & A, const float3 & x)
{
	return make_float3(A.m00*x.x + A.m01*x.y + A.m02*x.z,
					   A.m10*x.x + A.m11*x.y + A.m12*x.z,
					   A.m20*x.x + A.m21*x.y + A.m22*x.z);
}

__device__ float3x3 operator-(const float3x3 & a, const float3x3 & b)
{
	float3x3 c=a;
	for(int ii = 0;ii<float3x3::MATLEN;ii++){
		c.m[ii] -= b.m[ii];
	}
	return c;
}

__device__ float3x3 operator+(const float3x3 & a, const float3x3 & b)
{
	float3x3 c=a;
	for(int ii = 0;ii<float3x3::MATLEN;ii++){
		c.m[ii] += b.m[ii];
	}
	return c;
}
#endif