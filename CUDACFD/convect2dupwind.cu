/*
 * convect2dupwind.cu
 * 1-dimensional convection with time-independent velocty vector field using
 * CUDA C/C++ implementing "Upwind" interpolation
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160620
*/
#include "commonlib/gpu_anim.h"

#include "cuda.h"
#include "math.h" // CUDA C/C++ math.h

#define NthreadsX 800
#define NthreadsY 800
#define DELTAt 1./(10000000.)
#define RHO0 0.656
#define L_0X 1.0
#define L_0Y 1.0
#define DELTAx L_0X/((float) NthreadsX)
#define DELTAy L_0Y/((float) NthreadsY)

__global__ void convect( float *rho, float *ux, float *uy) {
	// map from threadIdx/blockIdx to x grid position
	int k_x = threadIdx.x + blockIdx.x * blockDim.x;
	int k_y = threadIdx.y + blockIdx.y * blockDim.y;
	int offset = k_x + k_y*blockDim.y*gridDim.y;
	
	int left     = offset - 1;
	int right    = offset + 1;
	int down     = offset + NthreadsX;  // on a bitmap, up down is "backwards"
	int up       = offset - NthreadsX;

	if (k_x == 0) {
		left++;
	}
	if (k_x==(NthreadsX-1)) {
		right--; }
	if (k_y == 0) { up += NthreadsX; }
	if (k_y == (NthreadsY - 1 )) { down -= NthreadsX ; }

/*
	if (ux[offset] > 0. ) {
		const float flux_right = DELTAy*rho[k_x][k_y]*ux[k_x][k_y] ; }
	else {
		const float flux_right = DELTAy*rho[right][k_y]*ux[k_x][k_y]; }
	if (ux[left] > 0. ) {
		const float flux_left = (-1.)*DELTAy*rho[left][k_y]*ux[left][k_y] ; }
	else {
		const float flux_left = (-1.)*DELTAy*rho[k_x][k_y]*ux[left][k_y]; }
	if (uy[offset] > 0. ) {
		const float flux_up = DELTAx*rho[k_x][k_y]*uy[k_x][k_y] ; }
	else {
		const float flux_up = DELTAx*rho[k_x][up]*uy[k_x][k_y]; }
	if (uy[down] > 0. ) {
		const float flux_down = (-1.)*DELTAx*rho[k_x][down]*uy[k_x][down] ; }
	else {
		const float flux_down = (-1.)*DELTAx*rho[k_x][k_y]*uy[k_x][down]; }
*/
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

//__global__ void float_to_color2d( uchar4* optr, const float outSrc[NthreadsX][NthreadsY]) {
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
	
	dim3 grids(NthreadsX/16,NthreadsY/16);
	dim3 threads(16,16);
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
	GPUAnimBitmap bitmap( NthreadsX, NthreadsY, &data );
	data.bitmap = &bitmap;
	data.totalTime = 0;
	data.frames = 0;
	// END of GPU animation setup END
	cudaEventCreate( &data.start );
	cudaEventCreate( &data.stop  );
	
	float rho[NthreadsX*NthreadsY];
	float ux[NthreadsX*NthreadsY];
	float uy[NthreadsX*NthreadsY];


	cudaMalloc((void**)&data.dev_rho, NthreadsX*NthreadsY*sizeof(float));
	cudaMalloc((void**)&data.dev_ux,NthreadsX*NthreadsY*sizeof(float));
	cudaMalloc((void**)&data.dev_uy,NthreadsX*NthreadsY*sizeof(float));

	// initial conditions
	for (int j=0; j<NthreadsY; ++j) {
		for (int i=0;i<NthreadsX; ++i) {
			ux[i+NthreadsX*j] = 10.0; // meters/second
			uy[i+NthreadsX*j] = 10.0; // meters/second
		}
	}
	
	for (int j=0; j<NthreadsY; ++j) {
		for (int i=0; i<NthreadsX; ++i) {
			rho[i+NthreadsX*j] = gaussian2d( ((float) i)*DELTAx,((float) j)*DELTAy,
									RHO0, 1./sqrt(0.000001),0.25,0.25);
		}
	}
	
	cudaMemcpy( data.dev_rho, rho, NthreadsX*NthreadsY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( data.dev_ux, ux, NthreadsX*NthreadsY*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( data.dev_uy, uy, NthreadsX*NthreadsY*sizeof(float), cudaMemcpyHostToDevice);
	
	
	bitmap.anim_and_exit((void (*)(uchar4*,void*,int))anim_gpu, NULL);
}
