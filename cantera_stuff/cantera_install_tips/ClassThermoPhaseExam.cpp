/*
 * ClassThermoPhaseExam.cpp
 * This is a simple, complete program that creates object representing gas mixture and prints its temperature
 * cf. http://www.cantera.org/docs/sphinx/html/cxx-guide/thermo.html#example-program
 * Computing Thermodynamic Properties
 * Class ThermoPhase
 *  
 * Compilation tip: I did this
 * g++ -std=c++11 ClassThermoPhaseExam.cpp -o ClassThermoPhaseExam -lcantera -l pthread
 * otherwise I obtain collect2: error: ld returned 1 exit status
 * */
#include "cantera/thermo.h"
#include <iostream>

int main(int argc, char** argv)
{
	Cantera::ThermoPhase* gas = Cantera::newPhase("h2o2.cti", "ohmech");
	std::cout << gas->temperature() << std::endl;
	return 0;
}
