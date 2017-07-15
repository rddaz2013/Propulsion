/* R1grid.h
 * R1 under discretization (discretize functor) to a grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160710
 */
#ifndef __R1GRID_H__
#define __R1GRID_H__

class Grid1d
{
	public :
		int L_x;
		float l_x;
		float h_x;
		float *rho;
		float *ux;
		
		Grid1d(int L, float l);
		float gridpt_to_space(int index_in);

		virtual ~Grid1d();

};

float gaussian1d( float A, float b, float c, float x) ;		 

#endif // _R1GRID_H__ 
