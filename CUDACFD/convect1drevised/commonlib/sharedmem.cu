/* sharedmem.cu
 * shared memory to a 1-dim. grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160710
 */
#include "sharedmem.h"

#include "errors.h" // HANDLE_ERROR()


namespace sharedmem {


__device__ Sharedmem :: Sharedmem(const int L_in, const int M_in, const int r_in) :
	L( L_in), M( M_in), r(r_in) { 
	N = (L + M - 1)/M;	
	S = M + 2*r;
	}
	
/*
__device__ void xder( float* f, Sharedmem sh_mem, StencilCoeff c_nus, float &result) {
	int k_x { threadIdx.x + blockIdx.x * blockDim.x };

	int L { sh_mem.L };
	int M { sh_mem.M };
	int S { sh_mem.S };
	
	int r { sh_mem.r };
	
	if (k_x >= L) { return; }

	__shared__ float sh_f[ M + 2*r ];

	// local s_x
	int s_x = threadIdx.x + r;
	
	sh_f[s_x] = f[k_x];
	
	// halo cells
	if (threadIdx.x < r) {
		sh_f[ s_x-r ] = f[ k_x -r ] ; 
		
		sh_f[s_x + M] = f[ k_x + M];
	}

	__syncthreads();

	int NUS { c_nus.NUS };

	__shared__ float stencil[ NUS ][2];
		
	for (int nu = 0; nu < NUS; ++nu ) {
		stencil[nu][0] = sh_f[ s_x - (nu +1)];
		stencil[nu][1] = sh_f[ s_x + (nu +1)];
	};
	
	__syncthreads();
	
	result = dev_dirder( stencil, c_nus);
}
*/

__device__ void xder1( float* f, Sharedmem sh_mem, float &result) {
	int k_x = threadIdx.x + blockIdx.x * blockDim.x ;

	constexpr int NUS = 1;

	const int L { sh_mem.L };
	const int M { sh_mem.M };
//	const int S { sh_mem.S };
	
	const int r { sh_mem.r };
	
	if (k_x >= L) { return; }

//	doesn't work for arrays, this initialization with non constant expressions __shared__ float sh_f[ M+2*r ];
	extern __shared__ float sh_f[ ];

	// local s_x
	int s_x = threadIdx.x + r;
	
	sh_f[s_x] = f[k_x];
	
	// halo cells
	if (threadIdx.x < r) {
		sh_f[ s_x-r ] = f[ k_x -r ] ; 
		
		sh_f[s_x + M] = f[ k_x + M];
	}

	__syncthreads();

	__shared__ float stencil[ NUS ][2];
	
		
	for (int nu = 0; nu < NUS; ++nu ) {
		stencil[nu][0] = sh_f[ s_x - (nu +1)];
		stencil[nu][1] = sh_f[ s_x + (nu +1)];
	};
	
	__syncthreads();
	
	result = dev_dirder1( stencil, dev_cnus);
}

__device__ void xder2( float* f, Sharedmem sh_mem, float &result) {
	int k_x = threadIdx.x + blockIdx.x * blockDim.x ;

	constexpr int NUS = 2;

	const int L { sh_mem.L };
	const int M { sh_mem.M };
//	const int S { sh_mem.S };
	
	const int r { sh_mem.r };
	
	if (k_x >= L) { return; }

	extern __shared__ float sh_f[ ];

	// local s_x
	int s_x = threadIdx.x + r;
	
	sh_f[s_x] = f[k_x];
	
	// halo cells
	if (threadIdx.x < r) {
		sh_f[ s_x-r ] = f[ k_x -r ] ; 
		sh_f[s_x + M] = f[ k_x + M];
	}

	__syncthreads();

	__shared__ float stencil[ NUS ][2];
	
	for (int nu = 0; nu < NUS; ++nu ) {
		stencil[nu][0] = sh_f[ s_x - (nu +1)];
		stencil[nu][1] = sh_f[ s_x + (nu +1)];
	};
	
	__syncthreads();
	
	result = dev_dirder2( stencil, dev_cnus);
}

__device__ void xder3( float* f, Sharedmem sh_mem, float &result) {
	int k_x = threadIdx.x + blockIdx.x * blockDim.x ;

	constexpr int NUS = 3;

	const int L { sh_mem.L };
	const int M { sh_mem.M };
//	const int S { sh_mem.S };
	
	const int r { sh_mem.r };
	
	if (k_x >= L) { return; }

	extern __shared__ float sh_f[ ];

	// local s_x
	int s_x = threadIdx.x + r;
	
	sh_f[s_x] = f[k_x];
	
	// halo cells
	if (threadIdx.x < r) {
		sh_f[ s_x-r ] = f[ k_x -r ] ; 
		sh_f[s_x + M] = f[ k_x + M];
	}

	__syncthreads();

	__shared__ float stencil[ NUS ][2];
	
	for (int nu = 0; nu < NUS; ++nu ) {
		stencil[nu][0] = sh_f[ s_x - (nu +1)];
		stencil[nu][1] = sh_f[ s_x + (nu +1)];
	};
	
	__syncthreads();
	
	result = dev_dirder3( stencil, dev_cnus);
}

__device__ void xder4( float* f, Sharedmem sh_mem, float &result) {
	int k_x = threadIdx.x + blockIdx.x * blockDim.x ;

	constexpr int NUS = 4;

	const int L { sh_mem.L };
	const int M { sh_mem.M };
//	const int S { sh_mem.S };
	
	const int r { sh_mem.r };
	
	if (k_x >= L) { return; }

	extern __shared__ float sh_f[ ];

	// local s_x
	int s_x = threadIdx.x + r;
	
	sh_f[s_x] = f[k_x];
	
	// halo cells
	if (threadIdx.x < r) {
		sh_f[ s_x-r ] = f[ k_x -r ] ; 
		sh_f[s_x + M] = f[ k_x + M];
	}

	__syncthreads();

	__shared__ float stencil[ NUS ][2];
	
	for (int nu = 0; nu < NUS; ++nu ) {
		stencil[nu][0] = sh_f[ s_x - (nu +1)];
		stencil[nu][1] = sh_f[ s_x + (nu +1)];
	};
	
	__syncthreads();
	
	result = dev_dirder4( stencil, dev_cnus);
}

}
