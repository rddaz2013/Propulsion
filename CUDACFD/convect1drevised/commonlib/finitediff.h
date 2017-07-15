/* finitediff.h
 * finite difference methods on a grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160710
 */
#ifndef __FINITEDIFF_H__
#define __FINITEDIFF_H__ 

#include "errors.h"

// I fixed the size of dev_c to be 4 because I don't anticipate that we'd need more accuracy 
extern __constant__ float dev_cnus[4];

void set1DerivativeParameters(float h_x );

void set2DerivativeParameters(float h_x );

void set3DerivativeParameters(float h_x );

void set4DerivativeParameters(float h_x );


__device__ float dev_dirder1(float stencil[1][2], float c_nus[4]);

__device__ float dev_dirder2(float stencil[2][2], float c_nus[4]);

__device__ float dev_dirder3(float stencil[3][2], float c_nus[4]);

__device__ float dev_dirder4(float stencil[4][2], float c_nus[4]);


#endif // __FINITEDIFF_H__
