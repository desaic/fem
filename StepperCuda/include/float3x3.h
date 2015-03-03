#ifndef FLOAT3X3_H
#define FLOAT3X3_H
#include "vector_types.h"
#include "cuda_runtime.h"

#define EPSILON 1e-7
#define SWAP(x,y) temp=x;x=y;y=temp

///@brief row-major order
class float3x3
{
public:
	static const int MATLEN=9;
	union{
		struct {float m00, m01, m02;
				float m10, m11, m12;
				float m20, m21, m22;
		};
		float m[MATLEN];
	};
	
	__device__ float3x3(){}
	
	__device__ float3x3& operator=(const float3x3 & aa)
	{
		for(int ii = 0;ii<MATLEN;ii++){
			m[ii] = aa.m[ii];
		}
		return *this;
	}
	
	__device__ float det2(float m00, float m01, 
					  float m10, float m11)const
	{
		return m00*m11-m01*m10;
	}

	__device__ void inverse(float3x3 & mInv)const
	{
		float cofactor00 =  det2( m11, m12, m21, m22 );
		float cofactor01 = -det2( m10, m12, m20, m22 );
		float cofactor02 =  det2( m10, m11, m20, m21 );

		float cofactor10 = -det2( m01, m02, m21, m22 );
		float cofactor11 =  det2( m00, m02, m20, m22 );
		float cofactor12 = -det2( m00, m01, m20, m21 );

		float cofactor20 =  det2( m01, m02, m11, m12 );
		float cofactor21 = -det2( m00, m02, m10, m12 );
		float cofactor22 =  det2( m00, m01, m10, m11 );

		float determinant = m00 * cofactor00 + m01 * cofactor01 + m02 * cofactor02;

		if( abs( determinant ) < EPSILON ){
			return;
		}
		else{
			float oneDet = 1.0f / determinant;
			mInv.m00 = cofactor00 * oneDet;
			mInv.m01 = cofactor10 * oneDet;
			mInv.m02 = cofactor20 * oneDet;

			mInv.m10 = cofactor01 * oneDet;
			mInv.m11 = cofactor11 * oneDet;
			mInv.m12 = cofactor21 * oneDet;

			mInv.m20 = cofactor02 * oneDet;
			mInv.m21 = cofactor12 * oneDet;
			mInv.m22 = cofactor22 * oneDet;

			return;
		}
	}

	__device__ float det()const{
	  return  m00 * ( m11 * m22 - m12 * m21 )
			- m01 * ( m10 * m22 - m12 * m20 )
			+ m02 * ( m10 * m21 - m11 * m20 );
	}


	__device__ void setIdentity()
	{
		setZero();
		for(int ii = 0;ii<3;ii++){
			m[ii*3+ii]=1;
		}
	}

	__device__ void setZero()
	{
		for(int ii = 0;ii<MATLEN;ii++){
			m[ii] =0;
		}
	}

	__device__ float3x3& operator += ( const float3x3 & x )
	{
		for(int ii = 0;ii<MATLEN;ii++){
			m[ii] += x.m[ii];
		}
		return *this;
	}

	__device__ void outerProd(const float3 & a, const float3 & b)
	{
		m00 = a.x*b.x;
		m01 = a.x*b.y;
		m02 = a.x*b.z;

		m10 = a.y*b.x;
		m11 = a.y*b.y;
		m12 = a.y*b.z;

		m20 = a.z*b.x;
		m21 = a.z*b.y;
		m22 = a.z*b.z;
	}

	__device__ float trace()const
	{
		return m00 + m11 + m22;
	}

	__device__ void transpose()
	{
		float temp;
		SWAP(m01,m10);
		SWAP(m02,m20);
		SWAP(m21,m12);
	}
};

#endif