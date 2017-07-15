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
		int L;
		int M;
		int r;
		int N;
		int S;
		
		__device__ Sharedmem(const int L_in, const int M_in, const int r_in );
		
};

//__device__ void xder( float* f, Sharedmem sh_mem, StencilCoeff c_nus, float &result);

__device__ void xder1( float* f, Sharedmem sh_mem, float &result);

__device__ void xder2( float* f, Sharedmem sh_mem, float &result);


__device__ void xder3( float* f, Sharedmem sh_mem, float &result);


__device__ void xder4( float* f, Sharedmem sh_mem, float &result);



}

#endif // __SHAREDMEM_H__

		

		
		
