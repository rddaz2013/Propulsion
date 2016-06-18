/*
 * convect1d.cu
 * 1-dimensional convection with time-independent velocty vector field using
 * CUDA C/C++ 
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160612
*/
#include "commonlib/gpu_anim.h"

#include "cuda.h"
#include "math.h" // CUDA C/C++ math.h

#define Nthreads 1600
#define DELTAt 1./(10000.)
#define RHO0 0.656
#define L_0 1.0
#define DELTAx L_0/((float) Nthreads)
#define DIMY 1600

__global__ void convect( float* rho, float *u) {
	// map from threadIdx/blockIdx to x grid position
	int k_x = threadIdx.x + blockIdx.x * blockDim.x;
	
	int left     = k_x - 1;
	int leftleft = k_x - 2;

	// boundary conditions
	if (k_x == 0) {
		left++;
		leftleft = leftleft + 2;
	}
	else if (k_x==1) {
		leftleft++;
	}
	
	float rholeftnew;
	rholeftnew = rho[left]-( rho[k_x]*u[k_x]-rho[leftleft]*u[leftleft])*(1./(2.*DELTAx))*DELTAt;
	
	rho[k_x] = (-1.)*rholeftnew + rho[k_x]+rho[left]   
				-( rho[k_x]*u[k_x]-rho[left]*u[left] )*(2.*DELTAt/DELTAx) ;
	rho[left] = rholeftnew ;
	__syncthreads();
}

float gaussian( float x, float A, float k, float x_0) 
{
	return A*exp(-k*(x-x_0)*(x-x_0)); 
}

// globals needed by the update routine
struct DataBlock {
	float         *dev_rho;
	float         *dev_u;
	GPUAnimBitmap *bitmap;
	cudaEvent_t   start, stop;
	float         totalTime;
	float         frames;
};

__global__ void float_to_1dimplot( uchar4* optr, const float* outSrc) {
	// map from threadIdx/BlockIdx to pixel position
	int x = threadIdx.x + blockIdx.x*blockDim.x ;
	
	int ivalue = ((int) (outSrc[x]*((float) DIMY) ));
	
	for (int j = 0; j < DIMY ; ++j ) {
		int offset = x + j*gridDim.x * blockDim.x;
		if (j < ivalue ){
			optr[offset].x = 0;
			optr[offset].y = 255;
			optr[offset].z = 0;
			optr[offset].w = 255;
		} else {
			optr[offset].x = 255;
			optr[offset].y = 0;
			optr[offset].z = 0;
			optr[offset].w = 255;
		}
	}
}

void anim_gpu(uchar4* outputBitmap, DataBlock *d, int ticks) {
	cudaEventRecord( d-> start,0 );
	
	/* change the 1000 time steps per frame manually */
	for (int i = 0; i < 5; ++i ) {
		convect<<<Nthreads/40,40>>>( d->dev_rho, d->dev_u);
	}
	
	float_to_1dimplot<<<Nthreads/40,40>>>( outputBitmap, d->dev_rho);
	
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
	GPUAnimBitmap bitmap( Nthreads, DIMY, &data );
	data.bitmap = &bitmap;
	data.totalTime = 0;
	data.frames = 0;
	// END of GPU animation setup END
	cudaEventCreate( &data.start );
	cudaEventCreate( &data.stop  );
	
	float rho[Nthreads];
	float u[Nthreads];
	
	cudaMalloc((void**)&data.dev_rho, Nthreads*sizeof(float));
	cudaMalloc((void**)&data.dev_u,Nthreads*sizeof(float));
	
	for (int j=0; j<Nthreads; ++j) {
		u[j] = 1.0; // meters/second
	}
	
	for (int j=0; j<Nthreads; ++j) {
		rho[j] = gaussian( ((float) j)*DELTAx,RHO0, 1./sqrt(0.0001),0.25);
	}
	
	cudaMemcpy( data.dev_rho, rho, Nthreads*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( data.dev_u, u, Nthreads*sizeof(float), cudaMemcpyHostToDevice);
	
	bitmap.anim_and_exit((void (*)(uchar4*,void*,int))anim_gpu, NULL);
}
