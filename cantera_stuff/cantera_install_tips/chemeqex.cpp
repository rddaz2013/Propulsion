/*
 * chemeqex.cpp
 * Chemical Equilibrium Example Program
 * 
 * cf. http://www.cantera.org/docs/sphinx/html/cxx-guide/equil-example.html
 * 
 * Compiling tip: I did this:
 * g++ -std=c++11 chemeqex.cpp -o chemeqex -lcantera -L/usr/local/lib64 -lsundials_cvodes -lsundials_ida -lsundials_nvecserial -L/usr/local/lib -l pthread
 * and it worked.  Then I started removing flags and found this was the "minimum" (number of flags) needed:
 * g++ -std=c++11 chemeqex.cpp -o chemeqex -lcantera -l pthread
 * 
 * */
#include <iostream>
#include <cantera/thermo.h>

using namespace Cantera;

void equil_demo()
{
	std::auto_ptr<ThermoPhase> gas(newPhase("h2o2.cti","ohmech"));
	gas->setState_TPX(1500.0, 2.0*OneAtm, "O2:1.0, H2:3.0, AR:1.0");
	gas->equilibrate("TP");
	std::cout << gas->report() << std::endl;
}

int main()
{
	try {
		equil_demo();
	} catch (CanteraError& err) {
		std::cout << err.what() << std::endl;
	}
}
