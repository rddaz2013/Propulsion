/* module.cpp
 * Examples of R-modules in C++11/14 style
 * Ernest Yeung  ernestyalumni@gmail.com
 * 20160630
 */
#include <iostream>
//#include <vector> // not needed
#include <array>
#include <chrono>

class complex {
	double re, im; // representation: two doubles
	public:
		complex(double r, double i) :re{r}, im{i} {}
		complex(double r) :re{r}, im{0} {}
		complex() :re{0}, im{0} {}
			// construct complex from two scalars
			// construct complex from one scalar
			// default complex: {0,0}

		double real() const { return re; }
		void real(double d) { re=d; }
		double imag() const { return im; }
		void imag(double d) { im=d; }
/*
		complex& operator+=(complex z) { re+=z.re , im+=z.im; return ∗this; }
			// add to re and im
			// and return the result
		complex& operator−=(complex z) { re−=z.re , im−=z.im; return ∗this; }
		complex& operator∗=(complex);
		complex& operator/=(complex);
			// defined out-of-class somewhere
			// defined out-of-class somewhere
*/
};


// class form
// \mathbb{R}^3 over \mathbb{R}, i.e. R^3 over RR
class RR_RR3 {
	public:
		float entries[3]; 

		RR_RR3() : entries{0.f,0.f,0.f} {}
		RR_RR3(float x1, float x2, float x3) : entries{ x1, x2, x3 } {}
		RR_RR3(float x1, float x2) : entries{ x1, x2, 0.f } {}
		RR_RR3(float x1) : entries{ x1, 0.f, 0.f } {}
		
		float xi(const int i) const { return entries[i]; };
		float* xis() { 
			float *xisptr = entries; 
			return xisptr; };

		RR_RR3& operator+=(RR_RR3 y) {
			for (auto i = 0; i<3; ++i) {
				entries[i]+= y.xi(i); }
				return *this; }

	//	~RR_RR3() { delete[] entry; }
		RR_RR3& operator*=(float scalar) {
//			for (auto i = 0 ; i < 3; ++i) {
			for (auto&& entry : entries) {
				entry*=scalar;
			}
			return *this;
		}

		RR_RR3& operator*=(float rotMat[][3]) {
			const float temp[3] { entries[0], entries[1], entries[2] };
//			std::cout << std::endl;
			for (auto i = 0; i<3; ++i) {
				float rowsum {0.f};
				for (auto j = 0; j<3; ++j) {
					// sanity check
//					std::cout << " " << rotMat[i][j] << " * " << temp[j] << " + ";
					rowsum += rotMat[i][j]*temp[j];

				}
				entries[i] = rowsum;
//				std::cout << std::endl;
			}
			return *this;
		}
};
RR_RR3 operator+(RR_RR3 x, RR_RR3 y) { return x+=y; }
RR_RR3 operator*(float scalar, RR_RR3 x) { return x*=scalar; }
RR_RR3 operator*(float rotMat[][3], RR_RR3 x) { return x*= rotMat; }

/*
 * An array of vectors as a (trivial) vector bundle on \mathbb{R}^1 \equiv RR
 * */

template <class T, size_t N_FLAT>
using VectorBundle = T[N_FLAT];

/*
 * Play around with rotation matrices
 * */
/*
struct so3{
	const float   theta_is[3];
	
	using aRotMat = std::array<std::array<float,3>,3> ;
//	float RotMat() { 
*/

class so3 {
//	using aRotMat = std::array<std::array<float,3>,3> ;

//	using aRotMat = float[3][3];
	public:
		const float theta_is[3]; 
	
//		aRotMat L_x { {0.f,0.f,0.f},{0.f,0.f,-1.f},{0.f,1.f,0.f}}; 

		float L_x[3][3] { {0.f,0.f,0.f},{0.f,0.f,-1.f},{0.f,1.f,0.f}}; 
		float L_y[3][3] { {0.f,0.f,1.f},{0.f,0.f,0.f},{-1.f,0.f,0.f}}; 
		float L_z[3][3] { {0.f,-1.f,0.f},{1.f,0.f,0.f},{0.f,0.f,0.f}}; 


		so3() : theta_is{0.f,0.f,0.f} {}
		so3(float x1, float x2, float x3) : theta_is{ x1, x2, x3 } {}
		so3(float x1, float x2) : theta_is{ x1, x2, 0.f } {}
		so3(float x1) : theta_is{ x1, 0.f, 0.f } {}		

		float thetai(const int i) const { return theta_is[i]; } ;
		const float* thetais() { 
			const float *thetaisptr = theta_is; 
			return thetaisptr; };	

		
//		void thetaL(aRotMat targetL) {
		void thetaL(float targetL[][3]) {
			for (auto i = 0; i<3; ++i) {
				for (auto j = 0; j < 3 ; ++j ) {
					targetL[i][j] = theta_is[0]*L_x[i][j] + theta_is[1]*L_y[i][j] 
									+ theta_is[2]*L_z[i][j];
				}
			}
		}



		float** thetaL2() {
//			float targetL[3][3];
			float** targetL = 0;
			targetL = new float*[3];

			for (auto i = 0; i<3; ++i) {
				targetL[i] = new float[3];
				for (auto j = 0; j < 3 ; ++j ) {
					targetL[i][j] = theta_is[0]*L_x[i][j] + theta_is[1]*L_y[i][j] 
									+ theta_is[2]*L_z[i][j];
				}
			}
			return targetL;
		}

};	

int main() {
	RR_RR3 testRR_RR30 {};
	RR_RR3 testRR_RR31 {2.35482f}; // full width at half maximum of gaussian
	RR_RR3 testRR_RR32 {2.35482f,4.29193f}; // full width at tenth max of gaussian
	RR_RR3 testRR_RR33 {2.35482f, 4.29193f,0.886f}; // typical speed of Maxwell-Boltzmann distribution
//	float testptr0 = testRR_RR30::xi(1);
//	std::cout << testRR_RR30 << std::endl;
//	std::cout << testRR_RR31 << std::endl;
//	std::cout << testRR_RR32 << std::endl;

	float testval = testRR_RR30.xi(0); 
	std::cout << testval << std::endl;
	std::cout << " full width at half maximum of gaussian " << testRR_RR31.xi(0) << std::endl;
	
	auto test3ptr = testRR_RR33.xis();
	std::cout << " typical speed of Maxwell-Boltzmann distribution " << test3ptr[2] << std::endl;
	
//	for (auto entry : test3ptr) {
//		std::cout << " entry in test3ptr " << *entry << " on " << entry << std::endl;
	std::cout << " size of test3ptr : " << sizeof(test3ptr) << std::endl;

	for (auto i = 0; i < 3 ; ++i) {
		std::cout << test3ptr[i] << std::endl; }
	
	// demonstrating abelian group addition
	
	RR_RR3 testadd0 = testRR_RR30 + testRR_RR33;
	test3ptr = testadd0.xis();
	
	for (auto i = 0; i < 3 ; ++i) {
		std::cout << test3ptr[i] << std::endl; }

	RR_RR3 testmult0 = 3.f * testRR_RR32;
	
	test3ptr = testmult0.xis();
	
	for (auto i = 0; i < 3 ; ++i) {
		std::cout << test3ptr[i] << std::endl; }

// template argument deduction FAILED
//	auto whatistest3ptr = std::begin(test3ptr);
//	for (auto jentry : test3ptr) {
//		std::cout << jentry << std::endl; }

	RR_RR3 VectorBundle[690000]; // above 690000 Segmentation fault (core dumped)

	so3 rotMateg { 0.886f, 1.5f, 1.085f };

	std::cout << " Rotation matrix example, theta_1 : " << rotMateg.thetai(0) << std::endl;
//	rotMateg.thetai(0);

	auto lalg_rotMat { rotMateg.thetaL2() };

	for (auto i = 0; i<3; ++i) {
		for (auto j = 0; j< 3; ++j) {
			std::cout << " " << lalg_rotMat[i][j] << " " ;		
		}
	std::cout << std::endl;
	}

	float targetL[3][3];
	rotMateg.thetaL(targetL);

	std::cout << " \n This is the result of using a method where we obtain the Lie algebra rotation matrix element by passing in a multidimensional array but treated as a pointer ref \n";

	for (auto i = 0; i<3; ++i) {
		for (auto j = 0; j< 3; ++j) {
			std::cout << " " << targetL[i][j] << " " ;		
		}
	std::cout << std::endl;
	}


	for (auto i = 0; i<3; ++i) {
		std::cout << " " << testRR_RR33.xi(i) << " " ;		
	}
	
	auto Rtimesvtest = targetL * testRR_RR33; // matrix R times test vector v
	
	std::cout << "\n This is the result of multiplication by a rotation matrix" << std::endl;;

	for (auto i = 0; i<3; ++i) {
		std::cout << " " << Rtimesvtest.xi(i) << " " ;		
	}
	
	// sanity check
	RR_RR3 testRR_RR4 {}; 
	testRR_RR4 = Rtimesvtest;
	for (auto i = 0; i<3; ++i) {
		std::cout << " " << testRR_RR4.xi(i) << " " ;		
	}
	std::cout << " \n Compared with testRR_RR33 " << std::endl;

	for (auto i = 0; i<3; ++i) {
		std::cout << " " << testRR_RR33.xi(i) << " " ;		
	}


	testRR_RR33 = Rtimesvtest;

	std::cout << " \n After reassignment " << " " << std::endl;

	for (auto i = 0; i<3; ++i) {
		std::cout << " " << testRR_RR33.xi(i) << " " ;		
	}
	
	
	std::cout << " \n End of sanity check " << std::endl;
	
	
	auto start = std::chrono::steady_clock::now();
	
	for (auto&& vec : VectorBundle) {
		vec = targetL * testRR_RR33;
		testRR_RR33 = vec;
	}

	auto end = std::chrono::steady_clock::now();
	
	// Store the time difference between start and end
	auto diff = end-start;
	
	std::cout << std::chrono::duration <double> (diff).count() << " sec" << std::endl;

	for (auto i = 0; i<3; ++i) {
		std::cout << (VectorBundle[32]).xi(i) << " " ;
	}
}
