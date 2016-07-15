/* finitediff.cu
 * finite difference methods on a grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160707
 */
#include "finitediff.h"
#include "errors.h"

__constant__ float dev_cnus[4];

void set1DerivativeParameters(float h_x )
{
	
	float cnus[4] { (1.f/2.f)*1.f / h_x,0.f,0.f,0.f  } ;
	
	HANDLE_ERROR(
		cudaMemcpyToSymbol( dev_cnus, cnus, sizeof(float)*4, 0, cudaMemcpyHostToDevice) ); // offset from start is 0
			
}

void set2DerivativeParameters(float h_x )
{
	
	float cnus[4] { 2.f/3.f*(1.f/h_x), -1.f/12.f*(1.f/h_x), 0.f, 0.f} ;
	
	HANDLE_ERROR(
		cudaMemcpyToSymbol( dev_cnus, cnus, sizeof(float)*4, 0, cudaMemcpyHostToDevice) ); // offset from start is 0
			
}

void set3DerivativeParameters(float h_x )
{
	
	float cnus[4] {  3.f / 4.f  * (1.f / h_x ) , -3.f / 20.f * (1.f / h_x ) ,
					1.f / 60.f * (1.f / h_x ) , 0.f} ;
	
	HANDLE_ERROR(
		cudaMemcpyToSymbol( dev_cnus, cnus, sizeof(float)*4, 0, cudaMemcpyHostToDevice) ); // offset from start is 0
			
}

void set4DerivativeParameters(float h_x )
{
	
	float cnus[4] { 4.f / 5.f  * (1.f / h_x ) , -1.f / 5.f * (1.f / h_x ) ,
					4.f / 105.f * (1.f / h_x ) , -1.f/280.f * (1.f/h_x)} ;
	
	HANDLE_ERROR(
		cudaMemcpyToSymbol( dev_cnus, cnus, sizeof(float)*4, 0, cudaMemcpyHostToDevice) ); // offset from start is 0
			
}

__device__ float dev_dirder1(float stencil[1][2], float c_nus[4]) {
	float tempvalue {0.f};

	tempvalue += c_nus[0]*( stencil[0][1] - stencil[0][0] );

	return tempvalue;
}

__device__ float dev_dirder2(float stencil[2][2], float c_nus[4]) {
	int NU {2};
	float tempvalue {0.f};

	for (int nu = 0; nu < NU; ++nu ) {
		tempvalue += c_nus[nu]*( stencil[nu][1] - stencil[nu][0] );
	}
	return tempvalue;
}

__device__ float dev_dirder3(float stencil[3][2], float c_nus[4]) {
	int NU {3};
	float tempvalue {0.f};
		
	for (int nu = 0; nu < NU; ++nu ) {
		tempvalue += c_nus[nu]*( stencil[nu][1] - stencil[nu][0] );
	}
	return tempvalue;
}

__device__ float dev_dirder4(float stencil[4][2], float c_nus[4]) {
	int NU {4};
	float tempvalue {0.f};

	for (int nu = 0; nu < NU; ++nu ) {
		tempvalue += c_nus[nu]*( stencil[nu][1] - stencil[nu][0] );
	}
	return tempvalue;
}

