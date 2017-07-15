/*
 * convect2dupwind_draft.cu
 * 2-dimensional convection with time-independent velocty vector field using
 * CUDA C/C++ implementing "Upwind" interpolation
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160620
*/
#include "commonlib/gpu_anim.h"

#include "cuda.h"
#include "math.h" // CUDA C/C++ math.h

#define N_x 40 // total number of blocks in the grid in x-direction
#define N_y 40 // total number of blocks in the grid in y-direction

#define M_x 20 // number of threads per block in x-direction
#define M_y 20 // number of threads per block in y-direction

#define DELTAt 1./(10000000.)
#define RHO0 0.656
#define L_0X 1.0
#define L_0Y 1.0
#define DELTAx L_0X/((float) N_x*M_x)
#define DELTAy L_0Y/((float) N_y*M_y)


__global__ void convect( float *rho, float *ux, float *uy) {
	// map from threadIdx/blockIdx to x grid position
	int k_x = threadIdx.x + blockIdx.x * blockDim.x;
	int k_y = threadIdx.y + blockIdx.y * blockDim.y;
	int offset = k_x + k_y*blockDim.x*gridDim.x;
	
	int left     = offset - 1;
	int right    = offset + 1;
	int down     = offset + N_x*M_x;  // on a bitmap, up down is "backwards"
	int up       = offset - N_x*M_x;

	if (k_x == 0) {
		left++;
	}
	if (k_x==(N_x*M_x-1)) {
		right--; }
	if (k_y == 0) { up += N_x*M_x; }
	if (k_y == (N_y*M_y - 1 )) { down -= N_x*M_x ; }

	float flux_right;
	if (ux[offset]>0.) { flux_right = DELTAy*rho[offset]*ux[offset] ; }
	else { flux_right = DELTAy*rho[right]*ux[offset] ; }
	float flux_left;
	if (ux[left] > 0.) { flux_left = (-1.)*DELTAy*rho[left]*ux[left]; }
	else { flux_left = (-1.)*DELTAy*rho[offset]*ux[left] ; }
	
	float flux_up;
	if (uy[up] > 0.) { flux_up = (-1.)*DELTAx*rho[up]*uy[up] ; }
	else { flux_up = (-1.)*DELTAx*rho[offset]*uy[up] ; }
	
	float flux_down;
	if (uy[offset] >0.) { flux_down = DELTAx*rho[offset]*uy[offset]; }
	else { flux_down = DELTAx*rho[down]*uy[offset]; }
		
	rho[offset] += (-1.)*(DELTAt/(DELTAx*DELTAy))*( flux_right+flux_left+flux_up+flux_down);
	
}
	
float gaussian2d( float x, float y, float A, float k, float x_0, float y_0) 
{
	return A*exp(-k*((x-x_0)*(x-x_0)+(y-y_0)*(y-y_0))); 
}

// globals needed by the update routine
struct DataBlock {
	float         *dev_rho;
	float         *dev_ux;
	float         *dev_uy;
	GPUAnimBitmap *bitmap;
	cudaEvent_t   start, stop;
	float         totalTime;
	float         frames;
};

__global__ void float_to_color2d( uchar4* optr, const float* outSrc) {
	// map from threadIdx/BlockIdx to pixel position
	int x = threadIdx.x + blockIdx.x*blockDim.x ;
	int y = threadIdx.y + blockIdx.y*blockDim.y ;
	int offset = x + y*blockDim.x*gridDim.x;
	
	float value = outSrc[offset];
	
	// Be aware of the "hard-coded" (numerical) constants for 
	// maximum and minimum scalar values that'll be assigned white and black, respectively
	if (value < 0.0001 ) { value = 0; }
	else if (value > 1.0 ) { value = 1.; } 

	/* convert to long rainbow RGB* */
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
	
	dim3 grids(N_x,N_y);
	dim3 threads(M_x,M_y);
	/* change the 1000 time steps per frame manually */
	for (int i = 0; i < 80; ++i ) {
		convect<<<grids,threads>>>( d->dev_rho, d->dev_ux, d->dev_uy);
	}
	
	float_to_color2d<<<grids,threads>>>( outputBitmap, d->dev_rho);
	
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
	DataBlock data;
	GPUAnimBitmap bitmap( N_x*M_x, N_y*M_y, &data );
	data.bitmap = &bitmap;
	data.totalTime = 0;
	data.frames = 0;
	// END of GPU animation setup END
	cudaEventCreate( &data.start );
	cudaEventCreate( &data.stop  );
	
	float rho[N_x*M_x*N_y*M_y];
	float ux[N_x*M_x*N_y*M_y];
	float uy[N_x*M_x*N_y*M_y];


	cudaMalloc((void**)&data.dev_rho, N_x*M_x*N_y*M_y*sizeof(float));
	cudaMalloc((void**)&data.dev_ux,N_x*M_x*N_y*M_y*sizeof(float));
	cudaMalloc((void**)&data.dev_uy,N_x*M_x*N_y*M_y*sizeof(float));

	// initial conditions
	for (int j=0; j<N_y*M_y; ++j) {
		for (int i=0;i<N_x*M_x; ++i) {
			ux[i+N_x*M_x*j] = 10.0; // meters/second
			uy[i+N_x*M_x*j] = 10.0; // meters/second
		}
	}
	
	for (int j=0; j<N_y*M_y; ++j) {
		for (int i=0; i<N_x*M_x; ++i) {
			rho[i+N_x*M_x*j] = gaussian2d( ((float) i)*DELTAx,((float) j)*DELTAy,
									RHO0, 1./sqrt(0.000001),0.25,0.25);
		}
	}
	
	cudaMemcpy( data.dev_rho, rho, N_x*M_x*N_y*M_y*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( data.dev_ux, ux, N_x*M_x*N_y*M_y*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( data.dev_uy, uy, N_x*M_x*N_y*M_y*sizeof(float), cudaMemcpyHostToDevice);
	
	
	bitmap.anim_and_exit((void (*)(uchar4*,void*,int))anim_gpu, NULL);
}
