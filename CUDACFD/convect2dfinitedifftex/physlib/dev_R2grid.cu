/* dev_R2grid.cu
 * R2 under discretization (discretize functor) to a grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160630
 */
#include "dev_R2grid.h"

__constant__ int dev_Ld[2];

__host__ dev_Grid2d::dev_Grid2d( dim3 Ld_in) : Ld(Ld_in)
{
	HANDLE_ERROR(
		cudaMalloc((void**)&this->dev_rho, this->NFLAT()*sizeof(float) ) );
	HANDLE_ERROR(
		cudaMalloc((void**)&this->dev_E, this->NFLAT()*sizeof(float) ) );
	HANDLE_ERROR(
		cudaMalloc((void**)&this->dev_u, this->NFLAT()*sizeof(float2) ) );
	HANDLE_ERROR(
		cudaMalloc((void**)&this->dev_p, this->NFLAT()*sizeof(float2) ) );
	
}

/*
__host__ dev_Grid3d::~dev_Grid3d() {
	HANDLE_ERROR( 
		cudaFree( this->dev_rho ) );
	HANDLE_ERROR(
		cudaFree( this->dev_E ) );
	HANDLE_ERROR(
		cudaFree( this->dev_u ) );
	
}
* */

__host__ int dev_Grid2d :: NFLAT() {
	return Ld.x*Ld.y;
}	



__device__ dev_block3d :: dev_block3d(unsigned int N_x, unsigned int N_y, unsigned int N_z  ) :
	N_is {N_x,N_y,N_z} 
{}


__device__ int dev_block3d :: flatten(
							int i_x, int i_y, int i_z) {
	return i_x+i_y*N_is[0]+i_z*N_is[0]*N_is[1];

}	
	
