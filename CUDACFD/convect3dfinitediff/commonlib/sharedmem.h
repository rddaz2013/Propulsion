/* sharedmem.h
 * shared memory to a 1-dim. grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160710
 */
#ifndef __SHAREDMEM_H__
#define __SHAREDMEM_H__

#include "finitediff.h"

namespace sharedmem {

class Sharedmem
{
	public:
		int3 L;
		int3 M;
		int r;
		int3 N;
		int3 S;
		
		__device__ Sharedmem(const int L_in[3], const int M_in[3], const int r_in );
		

		
};

__device__ int idxClip( int idx, int idxMax) ;
		
__device__ int flatten(int k_x, int k_y, int k_z, int L_x, int L_y, int L_z) ;

__device__ float dev_div1( float* rho, float3* u, Sharedmem sh_mem ) ;

__device__ float dev_div2( float* rho, float3* u, Sharedmem sh_mem ) ;

__device__ float dev_div3( float* rho, float3* u, Sharedmem sh_mem ) ;

__device__ float dev_div4( float* rho, float3* u, Sharedmem sh_mem ) ;



}

#endif // __SHAREDMEM_H__

		

		
		
