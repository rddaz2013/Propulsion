/*
 * convect1dupwind.cu
 * 1-dimensional convection with time-independent velocty vector field using
 * CUDA C/C++ on GPU implementing "Upwind" interpolation scheme	
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160618
*/

#include "./commonlib/errors.h"
#include "./commonlib/gpu_anim.h"
#include "./commonlib/finitediff.h" // set*DerivativeParameters (1,2,3,4), dev_dirder* (1,2,3,4)
#include "./commonlib/sharedmem.h"

#include "./physlib/R1grid.h"

#define L_0 1.f 
#define DELTAx L_0/((float) 400)

__global__ void convectupwind( float* mavg, float* u, const int Nthreads ) {
//__global__ void convect( float* mavg, float* u, const int Nthreads ) {
	float DELTAt { 0.1f };

	int k_x = threadIdx.x + blockIdx.x * blockDim.x ;
	
	int left     = k_x - 1;
	int right    = k_x + 1;
	
	// boundary conditions
	if (k_x == 0) {
		left++;
	}
	else if (k_x == (Nthreads - 1)) {
		right--;
	}	
	float flux_right;
	if (u[k_x] > 0. ) {
		flux_right = mavg[k_x]*u[k_x]; }
	else {
		flux_right = mavg[right]*u[k_x]; }
	float flux_left;
	if (u[left] > 0.) {
		flux_left = mavg[left]*u[left]; }
	else {
		flux_left = mavg[k_x]*u[left]; }
	mavg[k_x] += (-1.)*(DELTAt/ DELTAx)*( flux_right - flux_left);

}

// convectfinitediff_naiveshared works
/*
__global__ void convectfinitediff_naiveshared( float* rho, float* rho_in, float* u, const int L_in) {
	float DELTAt { 0.000001f };
	__syncthreads();

	int k_x = threadIdx.x + blockIdx.x*blockDim.x;

	rho_in[k_x] = rho[k_x];
	
	__shared__ float teststencil[1][2];
	
	int left  = k_x - 1;
	int right = k_x + 1;
	
	// check boundary conditions
	if (k_x == 0 ) {
		left++; }
	if (k_x >= (L_in - 1) ) {
		right--;
	}
	teststencil[0][0] = rho_in[left]  * u[left];
	teststencil[0][1] = rho_in[right] * u[right];
	 
	
	float result { dev_dirder1(teststencil, dev_cnus ) };

	__syncthreads();

	
	rho[k_x] += DELTAt*(-1.) * result;
	__syncthreads();
}
*/

//__global__ void convect( float* rho, float* rho_in, float* u, const int L_in) {
__global__ void convect( float* rho, float* rho_in, float* u, const int L_in, const int M_in) {
	// sharedmem tests
	sharedmem::Sharedmem testsh_mem(L_in, M_in, 1);

	float DELTAt { 0.000001f };
	__syncthreads();

	int k_x = threadIdx.x + blockIdx.x*blockDim.x;

	rho_in[k_x] = rho[k_x];
	
	float result;
	
	__syncthreads();
	
	sharedmem::xder1( rho_in, testsh_mem, result) ;  
	
/*	
	__shared__ float teststencil[1][2];
	
	int left  = k_x - 1;
	int right = k_x + 1;
	
	// check boundary conditions
	if (k_x == 0 ) {
		left++; }
	if (k_x >= (L_in - 1) ) {
		right--;
	}
	teststencil[0][0] = rho_in[left]  * u[left];
	teststencil[0][1] = rho_in[right] * u[right];
	 
	
	float result { dev_dirder1(teststencil, dev_cnus ) };
*/
	__syncthreads();

	
	rho[k_x] += DELTAt*(-1.) * (result);
	__syncthreads();
}


	

// this works for sure for convection with finite difference
__global__ void convectfinitediff_global( float* rho, float* rho_in, float* u, const int L_in) {
	float DELTAt { 0.00001f };
	__syncthreads();
		
	int k_x = threadIdx.x + blockIdx.x*blockDim.x;

	rho_in[k_x] = rho[k_x];
	
	int left  = k_x - 1;
	int right = k_x + 1;
	
	// check boundary conditions
	if (k_x == 0 ) {
		left++; }
	if (k_x >= (L_in - 1) ) {
		right--;
	}

	// update step 
	__syncthreads();
	
	rho[k_x] += DELTAt*(-1.f/2.f)*1.f/ (DELTAx)*( rho_in[right]*u[right] - rho_in[left]*u[left] );

	__syncthreads();

}

struct DataBlock {
	float 				*dev_rho;
	float               *dev_rho_in;
	float 				*dev_u;
	GPUAnimBitmap 		*bitmap;
// The following are for time-keeping
	cudaEvent_t 		start, stop; 
	float 				totalTime;
	float 				frames;
};


__global__ void float_to_1dimplot( uchar4* optr, const float* outSrc, const int DIMY) {
	// map from threadIdx/BlockIdx to pixel position
	int x = threadIdx.x + blockIdx.x * blockDim.x ;
	
	int ivalue = ((int) (outSrc[x]*((float) DIMY) ));
	
	// remember optr is a pointer to the buffer that OpenGL and CUDA SHARES
	for (int j = 0 ; j < DIMY ; ++j ) {
		int offset = x + j *gridDim.x * blockDim.x ;
		if (j < ivalue ) {
			optr[offset].x = 0 ;
			optr[offset].y = 255 ;
			optr[offset].z = 0;
			optr[offset].w = 255;
		} else {
			optr[offset].x = 255 ;
			optr[offset].y = 0;
			optr[offset].z = 0;
			optr[offset].w = 255;
		}
	}
}


void anim_gpu(uchar4* outputBitmap, DataBlock *d, int ticks) {
	// set these manually
	const int DIMY { 400 };
	const int L_in { 1000 };
	const int M_in { 2 };
	const int N_x { (L_in + M_in -1 ) / M_in } ;

	cudaEventRecord( d->start, 0 ) ;
	
	/* change the 1000 time steps per frame manually */
	for (int i =0; i< 1000; ++i) {
//		convectupwind<<<N_x,M_in>>>( d->dev_rho, d->dev_u, L_in);
		convect<<<N_x,M_in>>>( d->dev_rho, d->dev_rho_in,d->dev_u, L_in, M_in);

	}
		
	float_to_1dimplot<<<N_x,M_in>>>(outputBitmap,d->dev_rho, DIMY);
	
	// Recording time for rough benchmarking only	
	cudaEventRecord( d->stop, 0) ;
	cudaEventSynchronize( d->stop ) ;
	float elapsedTime;
	cudaEventElapsedTime( &elapsedTime, d->start, d->stop);
	
	d->totalTime += elapsedTime;
	++d->frames;
	printf( "Average Time per frame:  %3.1f ms\n", d->totalTime/d->frames );
// END of Recording time for rough benchmarking only, END
}

int main() {
	

	const float RHO0 { 0.656f };
	const int L_x { 1000 };
	const float l_x { 1.f }; 
	Grid1d hostgrid1d( L_x, l_x );

	const int HEIGHT { 400 };
	
	const float h_x { hostgrid1d.h_x };
	
	// sanity check 
	std::cout << " This is the h_x grid step value: " << h_x << std::endl;
	std::cout << " This is the previous h_x grid step value : " << DELTAx << std::endl;




	DataBlock data;
	GPUAnimBitmap bitmap( L_x, HEIGHT, &data);
	data.bitmap = &bitmap;
	data.totalTime = 0;
	data.frames = 0;
	// END of GPU animation setup END
	cudaEventCreate( &data.start );
	cudaEventCreate( &data.stop  );
	
	
	cudaMalloc((void**)&data.dev_rho, hostgrid1d.L_x*sizeof(float));
	cudaMalloc((void**)&data.dev_rho_in, hostgrid1d.L_x*sizeof(float));

	cudaMalloc((void**)&data.dev_u, hostgrid1d.L_x*sizeof(float));

	set1DerivativeParameters( hostgrid1d.h_x);
	
	// fill in host memory for time-independent u 
	for(int j = 0; j< hostgrid1d.L_x ; ++j) {
		hostgrid1d.ux[j] = 1.0; // m/s
	}

	for (int j = 0; j < hostgrid1d.L_x; ++j) {
		hostgrid1d.rho[j] =  gaussian1d(RHO0, 0.25*hostgrid1d.l_x ,0.05, ((float) j)*hostgrid1d.h_x); 
	}

	cudaMemcpy( data.dev_rho, hostgrid1d.rho, hostgrid1d.L_x*sizeof(float), cudaMemcpyHostToDevice);	
	cudaMemcpy( data.dev_rho_in, hostgrid1d.rho, hostgrid1d.L_x*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy( data.dev_u, hostgrid1d.ux, hostgrid1d.L_x*sizeof(float), cudaMemcpyHostToDevice);
	
	
	bitmap.anim_and_exit( (void (*)(uchar4* , void*,int))anim_gpu, NULL);
}
