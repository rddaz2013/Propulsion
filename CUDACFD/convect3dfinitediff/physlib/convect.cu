/* convect.cu
 * 3-dim. convection by finite difference with shared memory
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160720
 */
#include "convect.h"

__constant__ float dev_Deltat[1] ;

// convect_fd_naive_sh - convection with finite difference and naive shared memory scheme
__global__ void convect_fd_naive_sh( float* dev_rho, float3* dev_u ) {
	const int NU = 2;
	
	// map from threadIdx/blockIdx to x grid position
	int k_x = threadIdx.x + blockIdx.x * blockDim.x;
	int k_y = threadIdx.y + blockIdx.y * blockDim.y;
	int k_z = threadIdx.z + blockIdx.z * blockDim.z;
	
	int k = k_x + k_y*blockDim.x*gridDim.x + k_z*blockDim.x*gridDim.x*blockDim.y*gridDim.y ;

	int3* stencilindicesplus  = new int3[ NU ] ;
	int3* stencilindicesminus = new int3[ NU ] ;

	for (int nu = 0; nu < NU; ++nu ) {
		stencilindicesplus[  nu ].x = k + (nu + 1) ; 
		stencilindicesminus[ nu ].x = k - (nu + 1) ; 
		stencilindicesplus[  nu ].y = k + (nu + 1)*dev_Ld[0] ; 
		stencilindicesminus[ nu ].y = k - (nu + 1)*dev_Ld[0] ; 
		stencilindicesplus[  nu ].z = k + (nu + 1)*dev_Ld[0]*dev_Ld[1] ; 
		stencilindicesminus[ nu ].z = k - (nu + 1)*dev_Ld[0]*dev_Ld[1] ; 
	}

	int XI = 0;

	// check boundary conditions
	for (int nu = 0; nu < NU; ++nu) {
		if (k_x == nu ) {
			XI = NU - nu ;
			for (int xi = 0; xi < XI; ++xi ) {
				stencilindicesminus[ NU - 1 - xi ].x += XI- xi ;  
			}
		}
	
		if (k_y == nu ) {
			XI = NU - nu ;
			for (int xi = 0; xi < XI; ++xi) { 
				stencilindicesminus[ NU - 1 - xi ].y += (XI - xi)*dev_Ld[0] ;
			}
		}
		
		if (k_z == nu) {
			XI = NU - nu ;
			for (int xi = 0 ; xi < XI; ++xi ) {
				stencilindicesminus[ NU - 1 - xi ].y += (XI - xi)*dev_Ld[0]*dev_Ld[1] ;  
			}
		}
	
		if (k_x == (dev_Ld[0] - (nu + 1) ) ) {
			XI = NU - nu ;
			for (int xi = 0; xi < XI; ++xi ) {
				stencilindicesplus[ NU - 1 - xi].x -= XI-xi ;
			}
		}
		
		if (k_y == (dev_Ld[1] - (nu + 1) ) ) {
			XI = NU - nu ;
			for (int xi = 0; xi < XI; ++xi ) {
				stencilindicesplus[ NU - 1 - xi].y -= (XI-xi)*dev_Ld[0] ;
			}
		}
		
		if (k_z == (dev_Ld[2] - (nu + 1) ) ) {
			XI = NU - nu ;
			for (int xi = 0; xi < XI; ++xi ) {
				stencilindicesplus[ NU - 1 - xi].z -= (XI-xi)*dev_Ld[0]*dev_Ld[1] ;
			}
		}
		
	}
	
	__shared__ float3 stencil[NU][2] ; 
	
	for (int nu = 0 ; nu < NU; ++nu ) {
		stencil[nu][0].x = dev_rho[stencilindicesminus[nu].x]*dev_u[stencilindicesminus[nu].x].x  ;
		stencil[nu][1].x = dev_rho[stencilindicesplus[nu].x]*dev_u[stencilindicesplus[nu].x].x  ;
		stencil[nu][0].y = dev_rho[stencilindicesminus[nu].y]*dev_u[stencilindicesminus[nu].y].y  ;
		stencil[nu][1].y = dev_rho[stencilindicesplus[nu].y]*dev_u[stencilindicesplus[nu].y].y  ;
		stencil[nu][0].z = dev_rho[stencilindicesminus[nu].z]*dev_u[stencilindicesminus[nu].z].z  ;
		stencil[nu][1].z = dev_rho[stencilindicesplus[nu].z]*dev_u[stencilindicesplus[nu].z].z  ;
	}	
	
	float div_value { dev_div2( stencil ) } ;
	
	__syncthreads();
	
	dev_rho[k] +=  dev_Deltat[0] * (-1.f) * div_value ;		
			
	__syncthreads();		
			
}


__global__ void convectfd_naiveshared2( float* rho, float3* u ) {
	float Deltat { 0.00001f };
	
	// map from threadIdx/blockIdx to x grid position
	int k_x = threadIdx.x + blockIdx.x * blockDim.x;
	int k_y = threadIdx.y + blockIdx.y * blockDim.y;
	int k_z = threadIdx.z + blockIdx.z * blockDim.z;

	int offset = k_x + k_y*blockDim.x*gridDim.x + k_z*blockDim.x*gridDim.x*blockDim.y*gridDim.y;

	int left      = offset - 1;
	int left2     = offset - 2;
	int right    = offset + 1;
	int right2    = offset + 2;

	int up     = offset + dev_Ld[0];  
	int up2     = offset + 2*dev_Ld[0];  

	int down       = offset - dev_Ld[0];
	int down2       = offset - 2*dev_Ld[0];

	int top      = offset + dev_Ld[0]*dev_Ld[1];
	int top2      = offset + 2*dev_Ld[0]*dev_Ld[1];

	int bottom   = offset - dev_Ld[0]*dev_Ld[1];
	int bottom2   = offset - 2*dev_Ld[0]*dev_Ld[1];


	if (k_x == 0) {
		left++;
		left2 += 2;
	}
	if (k_x == 1) {
		left2++;
	}
	if (k_x==(dev_Ld[0]-1)) {
		right--; 
		right2 -= 2 ; }
	if (k_x==(dev_Ld[0]-2)) {
		right2--; }

	if (k_y == 0) { down += dev_Ld[0];
		down2 += 2*dev_Ld[0];
		 }
	if (k_y == 1) { down2 += dev_Ld[0]; }

	if (k_y == (dev_Ld[1] - 1 )) { up -= dev_Ld[0] ; 
		up2 -= 2*dev_Ld[0] ;
		}
	if (k_y == (dev_Ld[1]- 2) ) { up2 -= dev_Ld[0] ; }


	if (k_z == 0) { 
		bottom += dev_Ld[0]*dev_Ld[1]; 
		bottom2 += 2*dev_Ld[0]*dev_Ld[1] ; 
		}
	if (k_z == 1) {
		bottom2 += dev_Ld[0]*dev_Ld[1] ; 
	}

	if (k_z == (dev_Ld[2] - 1) ) { 
		top -= dev_Ld[0]*dev_Ld[1] ; 
		top2 -= 2* dev_Ld[0]*dev_Ld[1] ;
		}
	if (k_z == (dev_Ld[2] -2 )  ){
		top2 -= dev_Ld[0]*dev_Ld[1] ;
	}



	__shared__ float3 stencil[2][2];

	stencil[0][0].x = rho[left] * u[left].x     ;
	stencil[0][1].x = rho[right] * u[right].x   ;
	stencil[1][0].x = rho[left2] * u[left2].x     ;
	stencil[1][1].x = rho[right2] * u[right2].x   ;
	
	
	stencil[0][0].y = rho[down] * u[down].y     ;
	stencil[0][1].y = rho[up] * u[up].y         ;
	stencil[1][0].y = rho[down2] * u[down2].y     ;
	stencil[1][1].y = rho[up2] * u[up2].y         ;


	stencil[0][0].z = rho[bottom] * u[bottom].z ;
	stencil[0][1].z = rho[top] * u[top].z       ;
	stencil[1][0].z = rho[bottom2] * u[bottom2].z ;
	stencil[1][1].z = rho[top2] * u[top2].z       ;


	float div_value { dev_div2( stencil ) } ;
	
	 
	rho[offset] += Deltat * (-1.f) * div_value ; 


}		




