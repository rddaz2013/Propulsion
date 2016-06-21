/*
 * convect3dupwind.cu
 * 3-dimensional convection with time-independent velocty vector field using
 * CUDA C/C++ implementing "Upwind" interpolation
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160620
*/
#include "./commonlib/gpu_anim.h"
#include "./commonlib/errors.h"
#include "./commonlib/gpu_bitmap.h"

#include "cuda.h"
#include "math.h" // CUDA C/C++ math.h

#define N_x 120 // total number of blocks in the grid in x-direction
#define N_y 120 // total number of blocks in the grid in y-direction
#define N_z 8 // total number of blocks in the grid in z-direction

#define M_x 16 // number of threads per block in x-direction
#define M_y 16 // number of threads per block in y-direction
#define M_z 4 // number of threads per block in z-direction

#define DELTAt 1./(100000.)
#define RHO0 0.956 // 0.656
#define L_0X 1.0
#define L_0Y 1.0
#define L_0Z 1.0
#define DELTAx L_0X/((float) N_x*M_x)
#define DELTAy L_0Y/((float) N_y*M_y)
#define DELTAz L_0Z/((float) N_z*M_z)


__global__ void convect( float *rho, float *ux, float *uy, float *uz) {
	// map from threadIdx/blockIdx to x grid position
	int k_x = threadIdx.x + blockIdx.x * blockDim.x;
	int k_y = threadIdx.y + blockIdx.y * blockDim.y;
	int k_z = threadIdx.z + blockIdx.z * blockDim.z;
	int offset = k_x + k_y*blockDim.x*gridDim.x + k_z*blockDim.x*gridDim.x*blockDim.y*gridDim.y;

	int left     = offset - 1;
	int right    = offset + 1;
	int down     = offset + N_x*M_x;  // on a bitmap, up down is "backwards"
	int up       = offset - N_x*M_x;
	int top      = offset + N_x*M_x*N_y*M_y;
	int bottom   = offset - N_x*M_x*N_y*M_y;

	if (k_x == 0) {
		left++;
	}
	if (k_x==(N_x*M_x-1)) {
		right--; }
	if (k_y == 0) { up += N_x*M_x; }
	if (k_y == (N_y*M_y - 1 )) { down -= N_x*M_x ; }
	if (k_z == 0) { bottom += N_x*M_x*N_y*M_y; }
	if (k_z == (N_z*M_z - 1) ) { top -= N_x*M_x*N_y*M_y ; }

	float flux_right;
	if (ux[offset]>0.) { flux_right = DELTAy*DELTAz*rho[offset]*ux[offset] ; }
	else { flux_right = DELTAy*DELTAz*rho[right]*ux[offset] ; }
	float flux_left;
	if (ux[left] > 0.) { flux_left = (-1.)*DELTAy*DELTAz*rho[left]*ux[left]; }
	else { flux_left = (-1.)*DELTAy*DELTAz*rho[offset]*ux[left] ; }
	
	float flux_up;
	if (uy[up] > 0.) { flux_up = (-1.)*DELTAx*DELTAz*rho[up]*uy[up] ; }
	else { flux_up = (-1.)*DELTAx*DELTAz*rho[offset]*uy[up] ; }
	
	float flux_down;
	if (uy[offset] >0.) { flux_down = DELTAx*DELTAz*rho[offset]*uy[offset]; }
	else { flux_down = DELTAx*DELTAz*rho[down]*uy[offset]; }
		
	float flux_top;
	if (uz[offset]>0.) { flux_top = DELTAx*DELTAy*rho[offset]*uz[offset] ; }
	else { flux_top = DELTAx*DELTAy*rho[top]*uz[offset]; }
	float flux_bottom;
	if (uz[bottom]>0. ) { flux_bottom = (-1.)*DELTAx*DELTAy*rho[bottom]*uz[bottom] ; }
	else { flux_bottom = (-1.)*DELTAx*DELTAy*rho[offset]*uz[bottom] ; }
		
	rho[offset] += (-1.)*(DELTAt/(DELTAx*DELTAy*DELTAz))*( 
					flux_right+flux_left+flux_up+flux_down +flux_top+flux_bottom);

//	rho[offset] = rho[offset];
}
	
float gaussian3d( float x, float y, float z, float A, float k, float x_0, float y_0, float z_0) 
{
	return A*exp(-k*((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0) +(z-z_0)*(z-z_0))); 
}


// globals needed by the update routine
struct DataBlock {
	float         *dev_rho;
	float         *dev_ux;
	float         *dev_uy;
	float         *dev_uz;
	GPUAnimBitmap *bitmap;
	cudaEvent_t   start, stop;
	float         totalTime;
	float         frames;
};


// test DataBlock, Datablock2
struct DataBlock2 {
	float         *dev_rho;
	float         *dev_ux;
	float         *dev_uy;
	float         *dev_uz;
//	GPUAnimBitmap *bitmap;
	cudaEvent_t   start, stop;
	float         totalTime;
	float         frames;
};


__global__ void profile2d( uchar4* optr , const float* outSrc) {
	// map from threadIdx/BlockIdx to pixel position
	int k_x = threadIdx.x + blockIdx.x*blockDim.x ;
	int k_y = threadIdx.y + blockIdx.y*blockDim.y ;

	// choose at which z coordinate to make the slice in x-y plane
	int zcoordslice = blockDim.z*gridDim.z/2*1; 
	int offset = k_x + k_y*blockDim.x*gridDim.x ;
	int fulloffset = offset + zcoordslice*blockDim.x*gridDim.x*blockDim.y*gridDim.y ;
	float value = outSrc[fulloffset];
	
	if (value < 0.0001 ) { value = 0; }
	else if (value > 1.0 ) { value = 1.; } 

	// convert to long rainbow RGB* 
	value = value/0.20;
	int valueint  = ((int) floorf( value )); // this is the integer part 
	int valuefrac = ((int) floorf(255*(value-valueint)) );
	
	switch( valueint )
	{
		case 0:	optr[offset].x = 255; optr[offset].y = valuefrac; optr[offset].z = 0; 
		optr[offset].w = 255; 
		break;
		case 1:	optr[offset].x = 255-valuefrac; optr[offset].y = 255; optr[offset].z = 0; 
		optr[offset].w = 255; 
		break;
		case 2:	optr[offset].x = 0; optr[offset].y = 255; optr[offset].z = valuefrac; 
		optr[offset].w = 255; 
		break;
		case 3:	optr[offset].x = 0; optr[offset].y = 255-valuefrac; optr[offset].z = 255; 
		optr[offset].w = 255; 
		break;
		case 4:	optr[offset].x = valuefrac; optr[offset].y = 0; optr[offset].z = 255; 
		optr[offset].w = 255; 
		break;
		case 5:	optr[offset].x = 255; optr[offset].y = 0; optr[offset].z = 255; 
		optr[offset].w = 255; 
		break;
	}
}

__global__ void float_to_color3d( uchar4* optr, const float* outSrc) {
	// map from threadIdx/BlockIdx to pixel position
	int k_x = threadIdx.x + blockIdx.x*blockDim.x ;
	int k_y = threadIdx.y + blockIdx.y*blockDim.y ;

	// choose at which z coordinate to make the slice in x-y plane
	int zcoordslice = blockDim.z*gridDim.z/2*1; 

	int offset = k_x + k_y*blockDim.x*gridDim.x ;
	int fulloffset = offset + zcoordslice*blockDim.x*gridDim.x*blockDim.y*gridDim.y ;
	float value = outSrc[fulloffset];
	
	// Be aware of the "hard-coded" (numerical) constants for 
	// maximum and minimum scalar values that'll be assigned white and black, respectively
	if (value < 0.0001 ) { value = 0; }
	else if (value > 1.0 ) { value = 1.; } 

	// convert to long rainbow RGB* 
	value = value/0.20;
	int valueint  = ((int) floorf( value )); // this is the integer part 
	int valuefrac = ((int) floorf(255*(value-valueint)) );
	
	switch( valueint )
	{
		case 0:	optr[offset].x = 255; optr[offset].y = valuefrac; optr[offset].z = 0; 
		optr[offset].w = 255; 
		break;
		case 1:	optr[offset].x = 255-valuefrac; optr[offset].y = 255; optr[offset].z = 0; 
		optr[offset].w = 255; 
		break;
		case 2:	optr[offset].x = 0; optr[offset].y = 255; optr[offset].z = valuefrac; 
		optr[offset].w = 255; 
		break;
		case 3:	optr[offset].x = 0; optr[offset].y = 255-valuefrac; optr[offset].z = 255; 
		optr[offset].w = 255; 
		break;
		case 4:	optr[offset].x = valuefrac; optr[offset].y = 0; optr[offset].z = 255; 
		optr[offset].w = 255; 
		break;
		case 5:	optr[offset].x = 255; optr[offset].y = 0; optr[offset].z = 255; 
		optr[offset].w = 255; 
		break;
	}
}

void anim_gpu(uchar4* outputBitmap, DataBlock *d, int ticks) {
	cudaEventRecord( d-> start,0 );
	
	dim3 grids(N_x,N_y,N_z);
	dim3 threads(M_x,M_y,M_z);
//	/* change the 1000 time steps per frame manually 
	for (int i = 0; i < 2; ++i ) {
		convect<<<grids,threads>>>( d->dev_rho, d->dev_ux, d->dev_uy, d->dev_uz);
	}
	
	float_to_color3d<<<grids,threads>>>( outputBitmap, d->dev_rho);
	
	// Recording time for rough benchmarking, only
	cudaEventRecord( d-> stop, 0);
	cudaEventSynchronize(d-> stop);
	float elapsedTime;
	cudaEventElapsedTime( &elapsedTime, d->start, d->stop);
	
	d->totalTime += elapsedTime;
	++d->frames;
	printf("Average Time per frame:  %3.1f ms\n",d->totalTime/d->frames );
	// END of Recording time for rough benchmarking, only, END
}

int main( void ) {
//	DataBlock2 data;
//	GPUBitmap bitmap2( N_x*M_x, N_y*M_y );
	
	DataBlock data;

	GPUAnimBitmap bitmap( N_x*M_x, N_y*M_y, &data );
	data.bitmap = &bitmap;
	data.totalTime = 0;
	data.frames = 0;
	// END of GPU animation setup END

	cudaEventCreate( &data.start );
	cudaEventCreate( &data.stop  );
	
/*	float rho[N_x*M_x*N_y*M_y*N_z*M_z];
	float ux[N_x*M_x*N_y*M_y*N_z*M_z];
	float uy[N_x*M_x*N_y*M_y*N_z*M_z];
	float uz[N_x*M_x*N_y*M_y*N_z*M_z];
*/
	float *rho = (float *)malloc( sizeof(float)*N_x*M_x*N_y*M_y*N_z*M_z); 
	float *ux  = (float *)malloc( sizeof(float)*N_x*M_x*N_y*M_y*N_z*M_z); 
	float *uy  = (float *)malloc( sizeof(float)*N_x*M_x*N_y*M_y*N_z*M_z); 
	float *uz  = (float *)malloc( sizeof(float)*N_x*M_x*N_y*M_y*N_z*M_z); 


	HANDLE_ERROR( 
		cudaMalloc((void**)&data.dev_rho,N_x*M_x*N_y*M_y*N_z*M_z*sizeof(float))
		);
	HANDLE_ERROR(
		cudaMalloc((void**)&data.dev_ux,N_x*M_x*N_y*M_y*N_z*M_z*sizeof(float))
			);
	HANDLE_ERROR(
		cudaMalloc((void**)&data.dev_uy,N_x*M_x*N_y*M_y*N_z*M_z*sizeof(float))
		);
	HANDLE_ERROR(
		cudaMalloc((void**)&data.dev_uz,N_x*M_x*N_y*M_y*N_z*M_z*sizeof(float))
		);

	// initial conditions

	for (int k=0; k<(N_z*M_z); ++k) {
		for (int j=0; j<(N_y*M_y); ++j) {
			for (int i=0;i<(N_x*M_x); ++i) {
				ux[i+(N_x*M_x)*j+(N_x*M_x)*(N_y*M_y)*k] = 20.0; // meters/second
				uy[i+(N_x*M_x)*j+(N_x*M_x)*(N_y*M_y)*k] = 20.0; // meters/second
				uz[i+(N_x*M_x)*j+(N_x*M_x)*(N_y*M_y)*k] = 16.0 ;
			}
		}
	}


//	for (int k=0; k<NthreadsZ; ++k) {
//		for (int j=0; j<NthreadsY; ++j) {
//			for (int i=0; i<NthreadsX; ++i) {
//				rho[i+NthreadsX*j+NthreadsX*NthreadsY*k] = 0.5 ; 
//			}
//		}
//	}
	
	for (int k=0; k<(N_z*M_z); ++k) {
		for (int j=0; j<(N_y*M_y); ++j) {
			for (int i=0; i<(N_x*M_x); ++i) {
				rho[i+(N_x*M_x)*j+(N_x*M_x)*(N_y*M_y)*k] = 
					gaussian3d( ((float) i)*DELTAx,((float) j)*DELTAy, ((float) k)*DELTAz,
									RHO0, 1./sqrt(0.0001),0.25,0.25, 0.5);
			}
		}
	}
	
	HANDLE_ERROR(
		cudaMemcpy( data.dev_rho, rho, N_x*M_x*N_y*M_y*N_z*M_z*sizeof(float), cudaMemcpyHostToDevice)
		);
	HANDLE_ERROR(
		cudaMemcpy( data.dev_ux, ux, N_x*M_x*N_y*M_y*N_z*M_z*sizeof(float), cudaMemcpyHostToDevice)
		);
	HANDLE_ERROR(
		cudaMemcpy( data.dev_uy, uy, N_x*M_x*N_y*M_y*N_z*M_z*sizeof(float), cudaMemcpyHostToDevice)
		);
	HANDLE_ERROR(
		cudaMemcpy( data.dev_uz, uz, N_x*M_x*N_y*M_y*N_z*M_z*sizeof(float), cudaMemcpyHostToDevice)
		);
	
	free(rho);
	free(ux);
	free(uy);
	free(uz);
	
	
	bitmap.anim_and_exit((void (*)(uchar4*,void*,int))anim_gpu, NULL);
/*
 * 	
	*/
	
//	dim3 grids(N_x,N_y,N_z);
//	dim3 threads(M_x,M_y,M_z);
//	/* change the 1000 time steps per frame manually 
//	for (int i = 0; i < 4; ++i ) {
//		convect<<<grids,threads>>>( data.dev_rho, data.dev_ux, data.dev_uy, data.dev_uz);
//	}

	/*
	cudaFree(data.dev_rho);
	cudaFree(data.dev_ux);
	cudaFree(data.dev_uy);
	cudaFree(data.dev_uz);
//	*/
	
//	bitmap.free_resources();
//	return 0;
//	*/
	
	
//	profile2d<<<grids,threads>>>( bitmap2.devPtr , data.dev_rho);
//	bitmap2.display_and_exit(); 

}
