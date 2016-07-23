/* sharedmem.cu
 * shared memory to a 1-dim. grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160710
 */
#include "sharedmem.h"

#include "errors.h" // HANDLE_ERROR()


namespace sharedmem {


__device__ Sharedmem :: Sharedmem(const int L_in[3], const int M_in[3], const int r_in) :
	r(r_in) { 
	L.x = L_in[0] ; 
	L.y = L_in[1] ; 
	L.z = L_in[2] ; 
	
	M.x = M_in[0] ; 
	M.y = M_in[1] ; 
	M.z = M_in[2] ; 

	N.x = ( L.x + M.x - 1)/M.x;
	N.y = ( L.y + M.y - 1)/M.y;
	N.z = ( L.z + M.z - 1)/M.z;

	S.x = M.x + 2*r;
	S.y = M.y + 2*r;
	S.z = M.z + 2*r;
}
	
__device__ int idxClip( int idx, int idxMax) {
			return idx > (idxMax - 1) ? (idxMax - 1) : (idx < 0 ? 0 : idx);
}
		
__device__ int flatten(int k_x, int k_y, int k_z, int L_x, int L_y, int L_z) {
	return idxClip( k_x, L_x) + idxClip( k_y,L_y)*L_x + idxClip(k_z, L_z)*L_x*L_y ;
}	
	
__device__ void dev_div1( float* rho, float3* u, Sharedmem sh_mem, float &result ) {
	int k_x = threadIdx.x + blockIdx.x * blockDim.x; 
	int k_y = threadIdx.y + blockIdx.y * blockDim.y; 
	int k_z = threadIdx.z + blockIdx.z * blockDim.z; 
	
	constexpr int NUS = 1;

	const int3 L = sh_mem.L ;
//	const int3 M = sh_mem.M ;

	const int r { sh_mem.r };

	const int3 S { blockDim.x + 2*r , blockDim.y + 2*r , blockDim.z + 2*r } ;
//	S.x = blockDim.x + 2 * r ; 
//	S.y = blockDim.y + 2 * r ; 
//	S.z = blockDim.z + 2 * r ; 



// EY : 20160718 check this; I want nothing to be done, but I return 0.f	
	if (k_x >= L.x || k_y >= L.y || k_z >= L.z ) { return ;  }

	const int k = flatten( k_x,k_y,k_z,L.x,L.y,L.z);

	extern __shared__ float sh_rho[ ];
	extern __shared__ float3 sh_u[ ];


	// local s_i
	const int s_x = threadIdx.x + r;
	const int s_y = threadIdx.y + r;
	const int s_z = threadIdx.z + r;

	const int s_k = flatten( s_x,s_y,s_z, S.x, S.y, S.z );

	// Load regular cells
	sh_rho[s_k] = rho[k];
	sh_u[s_k]   = u[k];
	
	// Load halo cells
	if (threadIdx.x < r) {
		sh_rho[ flatten( s_x -r,s_y,s_z,S.x,S.y,S.z)] = rho[ flatten(k_x-r,k_y,k_z,L.x,L.y,L.z)] ;
		sh_u[ flatten( s_x -r,s_y,s_z,S.x,S.y,S.z)] = u[ flatten(k_x-r,k_y,k_z,L.x,L.y,L.z)] ;
		
		sh_rho[flatten(s_x+blockDim.x,s_y,s_z,S.x,S.y,S.z)] = 
			rho[flatten(k_x+blockDim.x,k_y,k_z,L.x,L.y,L.z)];
		sh_u[flatten(s_x+blockDim.x,s_y,s_z,S.x,S.y,S.z)]   = 
			u[flatten(k_x+blockDim.x,k_y,k_z,L.x,L.y,L.z)]  ;
	}
	
	if (threadIdx.y < r) {
		sh_rho[ flatten( s_x,s_y-r,s_z,S.x,S.y,S.z)] = rho[ flatten(k_x,k_y-r,k_z,L.x,L.y,L.z)] ;
		sh_u[ flatten( s_x ,s_y-r,s_z,S.x,S.y,S.z)] = u[ flatten(k_x,k_y-r,k_z,L.x,L.y,L.z)] ;
		
		sh_rho[flatten(s_x,s_y+blockDim.y,s_z,S.x,S.y,S.z)] 
				= rho[flatten(k_x,k_y+blockDim.y,k_z,L.x,L.y,L.z)];
		sh_u[flatten(s_x,s_y+blockDim.y,s_z,S.x,S.y,S.z)]   = 
				u[flatten(k_x,k_y+blockDim.y,k_z,L.x,L.y,L.z)]  ;
	}

	if (threadIdx.z < r) {
		sh_rho[ flatten( s_x ,s_y,s_z-r,S.x,S.y,S.z)] = rho[ flatten(k_x,k_y,k_z-r,L.x,L.y,L.z)] ;
		sh_u[ flatten( s_x,s_y,s_z-r,S.x,S.y,S.z)] = u[ flatten(k_x,k_y,k_z-r,L.x,L.y,L.z)] ;
		
		sh_rho[flatten(s_x,s_y,s_z+blockDim.z,S.x,S.y,S.z)] = 
			rho[flatten(k_x,k_y,k_z+blockDim.z,L.x,L.y,L.z)];
		sh_u[flatten(s_x,s_y,s_z+blockDim.z,S.x,S.y,S.z)]   = 
			u[flatten(k_x,k_y,k_z+blockDim.z,L.x,L.y,L.z)]  ;
	}

	__syncthreads();

	__shared__ float3 stencil[NUS][2];
	
	
	for (int nu = 0; nu < NUS; ++nu ) {
		stencil[nu][0].x = sh_rho[flatten(s_x-(nu+1),s_y,s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x-(nu+1),s_y,s_z,S.x,S.y,S.z)].x;
		stencil[nu][1].x = sh_rho[flatten(s_x+(nu+1),s_y,s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x+(nu+1),s_y,s_z,S.x,S.y,S.z)].x;
		stencil[nu][0].y = sh_rho[flatten(s_x,s_y-(nu+1),s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y-(nu+1),s_z,S.x,S.y,S.z)].y;
		stencil[nu][1].y = sh_rho[flatten(s_x,s_y+(nu+1),s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y+(nu+1),s_z,S.x,S.y,S.z)].y;
		stencil[nu][0].z = sh_rho[flatten(s_x,s_y,s_z-(nu+1),S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y,s_z-(nu+1),S.x,S.y,S.z)].z;
		stencil[nu][1].z = sh_rho[flatten(s_x,s_y,s_z+(nu+1),S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y,s_z+(nu+1),S.x,S.y,S.z)].z;
	}
	float div_value { dev_div1( stencil ) };
	
//	return div_value;
	result = div_value;
}

	
__device__ float dev_div2( float* rho, float3* u, Sharedmem sh_mem ) {
	int k_x = threadIdx.x + blockIdx.x * blockDim.x; 
	int k_y = threadIdx.y + blockIdx.y * blockDim.y; 
	int k_z = threadIdx.z + blockIdx.z * blockDim.z; 
	
	constexpr int NUS = 2;

	const int3 L = sh_mem.L ;
//	const int3 M = sh_mem.M ;
	// const int3 S = sh_mem.S ;

	const int r { sh_mem.r };
	
	const int3 S { blockDim.x + 2*r, blockDim.y + 2*r , blockDim.z +2*r };
//	S.x = blockDim.x + 2 * r ; 
//	S.y = blockDim.y + 2 * r ; 
//	S.z = blockDim.z + 2 * r ; 
	
	// EY : 20160718 check this; I want nothing to be done, but I return 0.f	
//	if (k_x >= L.x || k_y >= L.y || k_z >= L.z ) { return 0.f; }

	const int k = flatten( k_x,k_y,k_z,L.x,L.y,L.z);

	extern __shared__ float sh_rho[ ];
	extern __shared__ float3 sh_u[ ];


	// local s_i
	const int s_x = threadIdx.x + r;
	const int s_y = threadIdx.y + r;
	const int s_z = threadIdx.z + r;

	const int s_k = flatten( s_x,s_y,s_z, S.x, S.y, S.z );

	// Load regular cells
	sh_rho[s_k] = rho[k];
	sh_u[s_k]   = u[k];
	
	// Load halo cells
	if (threadIdx.x < r) {
		sh_rho[ flatten( s_x -r,s_y,s_z,S.x,S.y,S.z)] = rho[ flatten(k_x-r,k_y,k_z,L.x,L.y,L.z)] ;
		sh_u[ flatten( s_x -r,s_y,s_z,S.x,S.y,S.z)] = u[ flatten(k_x-r,k_y,k_z,L.x,L.y,L.z)] ;
		
//		sh_rho[flatten(s_x+M.x,s_y,s_z,S.x,S.y,S.z)] = rho[flatten(k_x+M.x,k_y,k_z,L.x,L.y,L.z)];
//		sh_u[flatten(s_x+M.x,s_y,s_z,S.x,S.y,S.z)]   = u[flatten(k_x+M.x,k_y,k_z,L.x,L.y,L.z)]  ;
	
		sh_rho[flatten(s_x+blockDim.x,s_y,s_z,S.x,S.y,S.z)] = 
			rho[flatten(k_x+blockDim.x,k_y,k_z,L.x,L.y,L.z)];
		sh_u[flatten(s_x+blockDim.x,s_y,s_z,S.x,S.y,S.z)]   = 
			u[flatten(k_x+blockDim.x,k_y,k_z,L.x,L.y,L.z)]  ;
	}

	
	
	if (threadIdx.y < r) {
		sh_rho[ flatten( s_x,s_y-r,s_z,S.x,S.y,S.z)] = rho[ flatten(k_x,k_y-r,k_z,L.x,L.y,L.z)] ;
		sh_u[ flatten( s_x ,s_y-r,s_z,S.x,S.y,S.z)] = u[ flatten(k_x,k_y-r,k_z,L.x,L.y,L.z)] ;
		
//		sh_rho[flatten(s_x,s_y+M.y,s_z,S.x,S.y,S.z)] = rho[flatten(k_x,k_y+M.y,k_z,L.x,L.y,L.z)];
//		sh_u[flatten(s_x,s_y+M.y,s_z,S.x,S.y,S.z)]   = u[flatten(k_x,k_y+M.y,k_z,L.x,L.y,L.z)]  ;

		sh_rho[flatten(s_x,s_y+blockDim.y,s_z,S.x,S.y,S.z)] 
				= rho[flatten(k_x,k_y+blockDim.y,k_z,L.x,L.y,L.z)];
		sh_u[flatten(s_x,s_y+blockDim.y,s_z,S.x,S.y,S.z)]   = 
				u[flatten(k_x,k_y+blockDim.y,k_z,L.x,L.y,L.z)]  ;

	}

	if (threadIdx.z < r) {
		sh_rho[ flatten( s_x ,s_y,s_z-r,S.x,S.y,S.z)] = rho[ flatten(k_x,k_y,k_z-r,L.x,L.y,L.z)] ;
		sh_u[ flatten( s_x,s_y,s_z-r,S.x,S.y,S.z)] = u[ flatten(k_x,k_y,k_z-r,L.x,L.y,L.z)] ;
		
//		sh_rho[flatten(s_x,s_y,s_z+M.z,S.x,S.y,S.z)] = rho[flatten(k_x,k_y,k_z+M.z,L.x,L.y,L.z)];
//		sh_u[flatten(s_x,s_y,s_z+M.z,S.x,S.y,S.z)]   = u[flatten(k_x,k_y,k_z+M.z,L.x,L.y,L.z)]  ;

		sh_rho[flatten(s_x,s_y,s_z+blockDim.z,S.x,S.y,S.z)] = 
			rho[flatten(k_x,k_y,k_z+blockDim.z,L.x,L.y,L.z)];
		sh_u[flatten(s_x,s_y,s_z+blockDim.z,S.x,S.y,S.z)]   = 
			u[flatten(k_x,k_y,k_z+blockDim.z,L.x,L.y,L.z)]  ;

	}

	__syncthreads();

	__shared__ float3 stencil[NUS][2];
	
	
	for (int nu = 0; nu < NUS; ++nu ) {
		stencil[nu][0].x = sh_rho[flatten(s_x-(nu+1),s_y,s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x-(nu+1),s_y,s_z,S.x,S.y,S.z)].x ;
		stencil[nu][1].x = sh_rho[flatten(s_x+(nu+1),s_y,s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x+(nu+1),s_y,s_z,S.x,S.y,S.z)].x ;
		stencil[nu][0].y = sh_rho[flatten(s_x,s_y-(nu+1),s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y-(nu+1),s_z,S.x,S.y,S.z)].y ;
		stencil[nu][1].y = sh_rho[flatten(s_x,s_y+(nu+1),s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y+(nu+1),s_z,S.x,S.y,S.z)].y ;
		stencil[nu][0].z = sh_rho[flatten(s_x,s_y,s_z-(nu+1),S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y,s_z-(nu+1),S.x,S.y,S.z)].z ;
		stencil[nu][1].z = sh_rho[flatten(s_x,s_y,s_z+(nu+1),S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y,s_z+(nu+1),S.x,S.y,S.z)].z ;
	}
	float div_value { dev_div2( stencil ) };
	
	return div_value;
}

	
__device__ float dev_div3( float* rho, float3* u, Sharedmem sh_mem ) {
	int k_x = threadIdx.x + blockIdx.x * blockDim.x; 
	int k_y = threadIdx.y + blockIdx.y * blockDim.y; 
	int k_z = threadIdx.z + blockIdx.z * blockDim.z; 
	
	constexpr int NUS = 3;

	const int3 L = sh_mem.L ;
	const int3 M = sh_mem.M ;
	const int3 S = sh_mem.S ;

	const int r { sh_mem.r };

// EY : 20160718 check this; I want nothing to be done, but I return 0.f	

	if (k_x >= L.x || k_y >= L.y || k_z >= L.z ) { return 0.f; }

	const int k = flatten( k_x,k_y,k_z,L.x,L.y,L.z);

	extern __shared__ float sh_rho[ ];
	extern __shared__ float3 sh_u[ ];


	// local s_i
	const int s_x = threadIdx.x + r;
	const int s_y = threadIdx.y + r;
	const int s_z = threadIdx.z + r;

	const int s_k = flatten( s_x,s_y,s_z, S.x, S.y, S.z );

	// Load regular cells
	sh_rho[s_k] = rho[k];
	sh_u[s_k]   = u[k];
	
	// Load halo cells
	if (threadIdx.x < r) {
		sh_rho[ flatten( s_x -r,s_y,s_z,S.x,S.y,S.z)] = rho[ flatten(k_x-r,k_y,k_z,L.x,L.y,L.z)] ;
		sh_u[ flatten( s_x -r,s_y,s_z,S.x,S.y,S.z)] = u[ flatten(k_x-r,k_y,k_z,L.x,L.y,L.z)] ;
		
		sh_rho[flatten(s_x+M.x,s_y,s_z,S.x,S.y,S.z)] = rho[flatten(k_x+M.x,k_y,k_z,L.x,L.y,L.z)];
		sh_u[flatten(s_x+M.x,s_y,s_z,S.x,S.y,S.z)]   = u[flatten(k_x+M.x,k_y,k_z,L.x,L.y,L.z)]  ;
	}
	
	if (threadIdx.y < r) {
		sh_rho[ flatten( s_x,s_y-r,s_z,S.x,S.y,S.z)] = rho[ flatten(k_x,k_y-r,k_z,L.x,L.y,L.z)] ;
		sh_u[ flatten( s_x ,s_y-r,s_z,S.x,S.y,S.z)] = u[ flatten(k_x,k_y-r,k_z,L.x,L.y,L.z)] ;
		
		sh_rho[flatten(s_x,s_y+M.y,s_z,S.x,S.y,S.z)] = rho[flatten(k_x,k_y+M.y,k_z,L.x,L.y,L.z)];
		sh_u[flatten(s_x,s_y+M.y,s_z,S.x,S.y,S.z)]   = u[flatten(k_x,k_y+M.y,k_z,L.x,L.y,L.z)]  ;
	}

	if (threadIdx.z < r) {
		sh_rho[ flatten( s_x ,s_y,s_z-r,S.x,S.y,S.z)] = rho[ flatten(k_x,k_y,k_z-r,L.x,L.y,L.z)] ;
		sh_u[ flatten( s_x,s_y,s_z-r,S.x,S.y,S.z)] = u[ flatten(k_x,k_y,k_z-r,L.x,L.y,L.z)] ;
		
		sh_rho[flatten(s_x,s_y,s_z+M.z,S.x,S.y,S.z)] = rho[flatten(k_x,k_y,k_z+M.z,L.x,L.y,L.z)];
		sh_u[flatten(s_x,s_y,s_z+M.z,S.x,S.y,S.z)]   = u[flatten(k_x,k_y,k_z+M.z,L.x,L.y,L.z)]  ;
	}

	__shared__ float3 stencil[NUS][2];
	
	
	for (int nu = 0; nu < NUS; ++nu ) {
		stencil[nu][0].x = sh_rho[flatten(s_x-(nu+1),s_y,s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x-(nu+1),s_y,s_z,S.x,S.y,S.z)].x ;
		stencil[nu][1].x = sh_rho[flatten(s_x+(nu+1),s_y,s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x+(nu+1),s_y,s_z,S.x,S.y,S.z)].x ;
		stencil[nu][0].y = sh_rho[flatten(s_x,s_y-(nu+1),s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y-(nu+1),s_z,S.x,S.y,S.z)].y ;
		stencil[nu][1].y = sh_rho[flatten(s_x,s_y+(nu+1),s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y+(nu+1),s_z,S.x,S.y,S.z)].y ;
		stencil[nu][0].z = sh_rho[flatten(s_x,s_y,s_z-(nu+1),S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y,s_z-(nu+1),S.x,S.y,S.z)].z ;
		stencil[nu][1].z = sh_rho[flatten(s_x,s_y,s_z+(nu+1),S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y,s_z+(nu+1),S.x,S.y,S.z)].z ;
	}
	float div_value { dev_div3( stencil ) };
	
	return div_value;
}

	
__device__ float dev_div4( float* rho, float3* u, Sharedmem sh_mem ) {
	int k_x = threadIdx.x + blockIdx.x * blockDim.x; 
	int k_y = threadIdx.y + blockIdx.y * blockDim.y; 
	int k_z = threadIdx.z + blockIdx.z * blockDim.z; 
	
	constexpr int NUS = 4;

	const int3 L = sh_mem.L ;
	const int3 M = sh_mem.M ;
	const int3 S = sh_mem.S ;

	const int r { sh_mem.r };

// EY : 20160718 check this; I want nothing to be done, but I return 0.f	
	
	if (k_x >= L.x || k_y >= L.y || k_z >= L.z ) { return 0.f; }

	const int k = flatten( k_x,k_y,k_z,L.x,L.y,L.z);

	extern __shared__ float sh_rho[ ];
	extern __shared__ float3 sh_u[ ];


	// local s_i
	const int s_x = threadIdx.x + r;
	const int s_y = threadIdx.y + r;
	const int s_z = threadIdx.z + r;

	const int s_k = flatten( s_x,s_y,s_z, S.x, S.y, S.z );

	// Load regular cells
	sh_rho[s_k] = rho[k];
	sh_u[s_k]   = u[k];
	
	// Load halo cells
	if (threadIdx.x < r) {
		sh_rho[ flatten( s_x -r,s_y,s_z,S.x,S.y,S.z)] = rho[ flatten(k_x-r,k_y,k_z,L.x,L.y,L.z)] ;
		sh_u[ flatten( s_x -r,s_y,s_z,S.x,S.y,S.z)] = u[ flatten(k_x-r,k_y,k_z,L.x,L.y,L.z)] ;
		
		sh_rho[flatten(s_x+M.x,s_y,s_z,S.x,S.y,S.z)] = rho[flatten(k_x+M.x,k_y,k_z,L.x,L.y,L.z)];
		sh_u[flatten(s_x+M.x,s_y,s_z,S.x,S.y,S.z)]   = u[flatten(k_x+M.x,k_y,k_z,L.x,L.y,L.z)] ;
	}
	
	if (threadIdx.y < r) {
		sh_rho[ flatten( s_x,s_y-r,s_z,S.x,S.y,S.z)] = rho[ flatten(k_x,k_y-r,k_z,L.x,L.y,L.z)] ;
		sh_u[ flatten( s_x ,s_y-r,s_z,S.x,S.y,S.z)] = u[ flatten(k_x,k_y-r,k_z,L.x,L.y,L.z)] ;
		
		sh_rho[flatten(s_x,s_y+M.y,s_z,S.x,S.y,S.z)] = rho[flatten(k_x,k_y+M.y,k_z,L.x,L.y,L.z)];
		sh_u[flatten(s_x,s_y+M.y,s_z,S.x,S.y,S.z)]   = u[flatten(k_x,k_y+M.y,k_z,L.x,L.y,L.z)] ;
	}

	if (threadIdx.z < r) {
		sh_rho[ flatten( s_x ,s_y,s_z-r,S.x,S.y,S.z)] = rho[ flatten(k_x,k_y,k_z-r,L.x,L.y,L.z)] ;
		sh_u[ flatten( s_x,s_y,s_z-r,S.x,S.y,S.z)] = u[ flatten(k_x,k_y,k_z-r,L.x,L.y,L.z)] ;
		
		sh_rho[flatten(s_x,s_y,s_z+M.z,S.x,S.y,S.z)] = rho[flatten(k_x,k_y,k_z+M.z,L.x,L.y,L.z)];
		sh_u[flatten(s_x,s_y,s_z+M.z,S.x,S.y,S.z)]   = u[flatten(k_x,k_y,k_z+M.z,L.x,L.y,L.z)] ;
	}

	__shared__ float3 stencil[NUS][2];
	
	
	for (int nu = 0; nu < NUS; ++nu ) {
		stencil[nu][0].x = sh_rho[flatten(s_x-(nu+1),s_y,s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x-(nu+1),s_y,s_z,S.x,S.y,S.z)].x ;
		stencil[nu][1].x = sh_rho[flatten(s_x+(nu+1),s_y,s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x+(nu+1),s_y,s_z,S.x,S.y,S.z)].x ;
		stencil[nu][0].y = sh_rho[flatten(s_x,s_y-(nu+1),s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y-(nu+1),s_z,S.x,S.y,S.z)].y ;
		stencil[nu][1].y = sh_rho[flatten(s_x,s_y+(nu+1),s_z,S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y+(nu+1),s_z,S.x,S.y,S.z)].y ;
		stencil[nu][0].z = sh_rho[flatten(s_x,s_y,s_z-(nu+1),S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y,s_z-(nu+1),S.x,S.y,S.z)].z ;
		stencil[nu][1].z = sh_rho[flatten(s_x,s_y,s_z+(nu+1),S.x,S.y,S.z)]*
							sh_u[flatten(s_x,s_y,s_z+(nu+1),S.x,S.y,S.z)].z ;
	}
	float div_value { dev_div4( stencil ) };
	
	return div_value;
}

}
