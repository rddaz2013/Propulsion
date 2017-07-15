/* R1grid.cpp
 * R1 under discretization (discretize functor) to a grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160710
 */
#include "R1grid.h"
#include <cstdlib> // exit
#include <cmath> // exp

Grid1d :: Grid1d(int L, float l) : 
L_x(L), l_x(l) 
{
	h_x = l_x/L_x;
	rho = new float[L_x];
	ux = new float[L_x];
}

float Grid1d::gridpt_to_space(int index_in) {
	if ( (index_in < 0 ) || index_in >= L_x ) {
		std::exit(0);
	}
	return ( index_in*h_x);
}

Grid1d::~Grid1d() {
	delete[] rho;
	delete[] ux;
}

// initial conditions routines

// gaussian1d - look up wikipedia, Normal distribution, for the formula
float gaussian1d( float A, float b, float c, float x) 
{
	return A*exp( - (x-b)*(x-b) / (2.f * c*c));
}


	

