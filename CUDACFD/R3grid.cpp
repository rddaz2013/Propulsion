/* R3grid.cpp
 * R3 under discretization (discretize functor) to a grid
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160630
 */
#include <iostream>
// #include <cstdlib> // not needed
#include <array>

struct spatialN{
	const float  l_is[3];
};

struct Grid3d{
	const int   N_is[3];
	const float h_is[3];
	float l_i(int i) { 
		return (*this).N_is[i] * (*this).h_is[i];	
	}
	int* is_grid_pt(int* I) {
		if (I[0] >= (*this).N_is[0] || I[1] >= (*this).N_is[1] || I[2] >= (*this).N_is[2]) {
//			std::exit(0); }
			return nullptr; }
//		else if (*i_y < (*this).N_is[1]) {
//			std::exit(0); }
		//		return const int[3] { *i_x, *i_y, *i_z }
		return I;
	}
//	float* grid_pt_to_space(int* I) {
	void gridpt_to_space(int* Iind, float* Xi) {
			auto checkingrid = (*this).is_grid_pt(Iind);
		if (checkingrid == nullptr) {
			std::exit(0); }
// sanity check
		std::cout << " x coordinate : " << Iind[0]*(*this).h_is[0] << std::endl;

		
		Xi[0] = Iind[0]*(*this).h_is[0]; 
		Xi[1] = Iind[1]*(*this).h_is[1]; 
		Xi[2] = Iind[2]*(*this).h_is[2];

//		return x_y_z;
	}
	// out to flatten
	long int N_xN_yN_z() {
		long int N {1};
		for (auto N_i : (*this).N_is) {
			N *= N_i; }
		return N;
	}
	
	int flatten(int i_x, int i_y, int i_z) {
		return i_x+i_y*(*this).N_is[0]+i_z*(*this).N_is[0]*(*this).N_is[1];
	}
};

/* C^{\infty}(\mathbb{R}^3) functions; C^k function space */
template <class T, size_t N_x, size_t N_y, size_t N_z>
using Ckfunc = std::array<std::array<std::array<T, N_z>,N_y>,N_x>; 

/*
struct CKfuncspace {
	Ckfunc& operator+= (Ckfunc) {
		return *this;	
	}
};
*/


template <class T, size_t N_FLAT>
using CkfuncspaceFLAT = T[N_FLAT];



// add to re and im
// and return the result


// eg means example; it's just an example 
// old school
const int   N_xeg = 200, N_yeg = 200, N_zeg = 50; 
// new school
const float h_xeg { 0.0001f }, h_yeg { 0.001f }, h_zeg { 0.01f }; 

 
int main() {
	const int   N_iseg[3] {N_xeg,N_yeg,N_zeg};  // N_is means N_i's
	const float h_iseg[3] {h_xeg,h_yeg,h_zeg};  // h_is means h_i's

	// new school range-over for
	for (auto N_i : N_iseg) {
		std::cout << " This is in N_i : " << N_i << std::endl;	
	}

// old school way of iterating over an array, treating arrays as a pointer
// cf. http://stackoverflow.com/questions/14595285/cannot-use-begin-or-end-on-an-array
	for (auto i = h_iseg ; i != h_iseg + 3; ++i)
		{
		std::cout << " This is in h_i : " << *i << " on " << i << std::endl;	
	}
	for (auto i = 0; i != 3; ++i) {
		std::cout << " This is i : " << i << " This is N_i, h_i, respectively " 
			<< N_iseg[i] << ' ' << h_iseg[i] << std::endl;
	}

	std::cout << " This is begin for N_i : " << *std::begin(N_iseg) << std::endl;
	std::cout << " This is end   for N_i : " << *std::end(N_iseg)   << std::endl;
/*	for (auto i = std::begin(h_is); i < std::end(h_is), ++i;) {
		std::cout << " This is i in h_i : " << i << std::endl;
	}*/
// cf. http://en.cppreference.com/w/cpp/language/range-for

	// cf. http://en.cppreference.com/w/cpp/iterator/next
	auto i_eg = std::begin(N_iseg);
	std::cout << " This begins N_i : " << *i_eg << std::endl;
	while (i_eg != std::end(N_iseg)) { 
		std::cout << " This is the next N_i : " << *i_eg << " on " << i_eg << std::endl;
		i_eg = std::next(i_eg,1);
	}
// cf. http://en.cppreference.com/w/cpp/iterator/advance
	auto j_eg = std::begin(h_iseg);
	std::cout << " This begins h_i : " << *j_eg << " on " << j_eg << std::endl;
	std::advance(j_eg,1);
	std::cout << " This is the last h_i : " << *j_eg << " on " << j_eg << std::endl;
// Notice how advance returns a void, UNLIKE next, which returns the "next" memory address (haha)
	j_eg = std::begin(h_iseg);
	while (j_eg != std::end(h_iseg)) {
		std::cout << " This is h_i, advanced by 1 : " << *j_eg << " on " << j_eg << std::endl;
		std::advance(j_eg,1);
	}
	// const_iterator cf. http://www.cplusplus.com/forum/beginner/126145/
	
	Grid3d grid3d = { {N_xeg, N_yeg, N_zeg}, {h_xeg, h_yeg, h_zeg} } ;
	
	std::cout << " This is grid3d's N's : " << grid3d.N_is[0] << 
		" " << grid3d.N_is[1] << " " << grid3d.N_is[2] << std::endl ;

	std::cout << " This is grid3d's h's : " << grid3d.h_is[0] << 
		" " << grid3d.h_is[1] << " " << grid3d.h_is[2] << std::endl ;

	std::cout << " This is grid3d's l's : " << grid3d.l_i(0) << 
		" " << grid3d.l_i(1) << " " << grid3d.l_i(2) << std::endl ;

	int p[3] {2,3,5} ;
	
	auto ptresult = grid3d.is_grid_pt(p);
	std::cout << 
		" This is the result of pt being in grid or not (1 being not a success; it's equal to a null pointer) : " 
			<< (ptresult == nullptr) << std::endl;

	float Xpt[3] {0.f,0.f,0.f};

	grid3d.gridpt_to_space(p,Xpt);
	
	std::cout << "Now Xpt is this: " << Xpt[0] << " " << Xpt[1] << " " << Xpt[2] << std::endl;

	std::cout << "\n First steps to implementing the flatten functor. " << std::endl;
	std::cout << "This is the (total) length of the flattened grid : " << grid3d.N_xN_yN_z() << std::endl;
	std::cout << "This is the index of the flattened point : " << grid3d.flatten(p[0],p[1],p[2]) << std::endl;

	Ckfunc<float, N_xeg, N_yeg, N_zeg> rho;

	constexpr int NFLATeg = N_xeg*N_yeg*N_zeg;
//	CkfuncspaceFLAT<float, 2000000> rhoFLAT; // Segmentation fault (core dumped)

//	float CkfuncspaceFLATnative[2000000]; // Segmentation fault (core dumped)

	float CkfuncFLATnative[88000]; // above 88000, Segmentation fault (core dumped)


}	

