/*
 * main.cu
 * 3-dimensional convection with time-independent velocty vector field using
 * CUDA C/C++ implementing finite difference
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160723
*/
#include <functional>

#include "./commonlib/bit_anim2d.h"
#include "./commonlib/errors.h"
#include "./commonlib/finitediff.h"
#include "./commonlib/sharedmem.h"

#include "./physlib/convect.h"
#include "./physlib/R3grid.h"
#include "./physlib/dev_R3grid.h"

#include "math.h" // CUDA C/C++ math.h

#define GL_GLEXT_PROTOTYPES // needed for identifier glGenBuffer, glBindBuffer, glBufferData, glDeleteBuffers

#include <GL/glut.h>

#include <cuda_runtime.h>

#include <cuda_gl_interop.h>


const float Deltat[1] { 0.00001f } ; 

// physics
const int W { 640 } ;
const int H { 640 } ;
const int DEPTH { 288 } ;

dim3 dev_L3 { static_cast<unsigned int>(W), 
				static_cast<unsigned int>(H),
				static_cast<unsigned int>(DEPTH) };
				
dev_Grid3d dev_grid3d( dev_L3 );		


__global__ void convect_finitediff( float* dev_rho, float3* dev_u ) {

	const int NU = 2;
	
	// map from threadIdx/blockIdx to x grid position
	int k_x = threadIdx.x + blockIdx.x * blockDim.x;
	int k_y = threadIdx.y + blockIdx.y * blockDim.y;
	int k_z = threadIdx.z + blockIdx.z * blockDim.z;

	int k = k_x + k_y*blockDim.x*gridDim.x + k_z*blockDim.x*gridDim.x*blockDim.y*gridDim.y ;

	__shared__ int3 stencilindicesplus[NU] ;
	__shared__ int3 stencilindicesminus[NU] ;

	for (int nu = 0; nu < NU; ++nu ) {
		stencilindicesplus[  nu ].x = k + (nu + 1) ; 
		stencilindicesminus[ nu ].x = k - (nu + 1) ; 
		stencilindicesplus[  nu ].y = k + (nu + 1)*dev_Ld[0] ; 
		stencilindicesminus[ nu ].y = k - (nu + 1)*dev_Ld[0] ; 
		stencilindicesplus[  nu ].z = k + (nu + 1)*dev_Ld[0]*dev_Ld[1] ; 
		stencilindicesminus[ nu ].z = k - (nu + 1)*dev_Ld[0]*dev_Ld[1] ; 
	}

	int XI = 0;

	for (int nu = 0; nu < NU; ++nu) {
		if (k_x == nu ) {
			XI = NU-nu;
			for (int xi = 0 ; xi < XI; ++xi ) {
				stencilindicesminus[ NU - 1 - xi].x += XI - xi ; 
			}
		}
		
		if (k_y == nu ) {
			XI = NU-nu;
			for (int xi = 0 ; xi < XI; ++xi ) {
				stencilindicesminus[ NU - 1 - xi].y += (XI - xi)*dev_Ld[0] ; 
			}
		}

		if (k_z == nu ) {
			XI = NU-nu;
			for (int xi = 0 ; xi < XI; ++xi ) {
				stencilindicesminus[ NU - 1 - xi].z += (XI - xi)*dev_Ld[0]*dev_Ld[1] ; 
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
	
	dev_rho[k] += dev_Deltat[0] * (-1.f) * div_value ; 

}

// graphics + physics

const dim3 M_i { 2, 2, 2 };

const int iters_per_render { 10 } ;

GPUAnim2dBit animbitmap( W, H );
GPUAnim2dBit* bitptr = &animbitmap;


void make_idle_func(int iters_per_render ) {
	int ticks = 1; 
	uchar4* devPtr; 
	size_t size;
	
	HANDLE_ERROR(
		cudaGraphicsMapResources(1, &(bitptr->cuda_pixbufferObj_resource), NULL ) );
	
	HANDLE_ERROR(
		cudaGraphicsResourceGetMappedPointer((void **)&devPtr, &size, 
			bitptr->cuda_pixbufferObj_resource ) 
	);
	
	ticks++ ;
	
	dim3 grids( (dev_L3.x+M_i.x-1)/M_i.x,(dev_L3.y+M_i.y-1)/M_i.y,(dev_L3.z+M_i.z-1)/M_i.z) ;
	
	for (int i = 0 ; i < iters_per_render; ++i ) {
		convect_finitediff<<<grids,M_i>>>(
			dev_grid3d.dev_rho, dev_grid3d.dev_u ); 
	}
	float_to_color3d<<<grids,M_i>>>(devPtr, dev_grid3d.dev_rho ) ;
	printf("Iteration complete : ticks %d \n ", ticks );
	

	HANDLE_ERROR(
		cudaGraphicsUnmapResources( 1, &(bitptr->cuda_pixbufferObj_resource), NULL ));
	glutPostRedisplay();
}; 

std::function<void()> idle_func = std::bind( make_idle_func, iters_per_render ) ;

void idle() {
	idle_func() ;
}

std::function<void()> display_func = std::bind( make_display_func, W,H) ;

void display() {
	display_func();
}


int main(int argc, char** argv) {
	// physics
	constexpr float rho_0 { 0.956 };
	constexpr std::array<int,3> LdS { W, H, DEPTH} ; 
	constexpr std::array<float,3> ldS { 1.f, 1.f, 1.f };

	HANDLE_ERROR(
		cudaMemcpyToSymbol( dev_Deltat, Deltat, sizeof(float)*1,0,cudaMemcpyHostToDevice) );

	const int Ld_to_const[3] { LdS[0], LdS[1], LdS[2] };
	
	HANDLE_ERROR(
		cudaMemcpyToSymbol( dev_Ld, Ld_to_const, sizeof(int)*3,0,cudaMemcpyHostToDevice) );


	Grid3d grid3d( LdS, ldS);

	const float hds[3] { grid3d.hd[0], grid3d.hd[1], grid3d.hd[2] };

	set2DerivativeParameters(hds );


	// graphics setup
	
	bitptr->initGLUT(&argc,argv) ;
	 
	glutKeyboardFunc( keyboard_func );
	glutMouseFunc( mouse_func );
	glutIdleFunc( idle );
	glutDisplayFunc( display) ;
	bitptr->initPixelBuffer();
		

	// initial conditions

	for (int k=0; k<(grid3d.Ld[2]); ++k) {
		for (int j=0; j<(grid3d.Ld[1]); ++j) {
			for (int i=0;i<(grid3d.Ld[0]); ++i) {

				grid3d.u[ grid3d.flatten(i,j,k) ].x = 25.0;  // meters/second
				grid3d.u[ grid3d.flatten(i,j,k) ].y = 25.0;  // meters/second
				grid3d.u[ grid3d.flatten(i,j,k) ].z = 12.0;  // meters/second

			}
		}
	}

	std::array<int,3> ix_in { 0, 0, 0 };
	std::array<float,3> b_0 { 0.25f*grid3d.ld[0], 0.25f*grid3d.ld[1], 0.5f*grid3d.ld[2]  };
	
	for (int k=0; k<(grid3d.Ld[2]); ++k) {
		for (int j=0; j<(grid3d.Ld[1]); ++j) {
			for (int i=0; i<(grid3d.Ld[0]); ++i) {
				ix_in[0] = i;
				ix_in[1] = j;
				ix_in[2] = k;
				grid3d.rho[ grid3d.flatten(i,j,k) ] = 
						gaussian3d( rho_0, 0.05, b_0,grid3d.gridpt_to_space(ix_in));

			}
		}
	}


	HANDLE_ERROR(
		cudaMemcpy( dev_grid3d.dev_rho, grid3d.rho, grid3d.NFLAT()*sizeof(float), cudaMemcpyHostToDevice)
		);

	HANDLE_ERROR(
		cudaMemcpy( dev_grid3d.dev_u, grid3d.u, grid3d.NFLAT()*sizeof(float3), cudaMemcpyHostToDevice)
		);


	glutMainLoop();

	HANDLE_ERROR(
		cudaFree( dev_grid3d.dev_rho )  );

	HANDLE_ERROR(
		cudaFree( dev_grid3d.dev_u ) );
		
	bitptr->exitfunc(); 
}
	
