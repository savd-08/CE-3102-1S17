#include <iostream>
#include <cstdlib>
#include "biseccion.h"
#include "interpolacion.h"


/**************************************************/

double check(const double t) {
	static const double pi = 3.14159265358979323856264338327950288;
	return 0.5*std::exp(-t) - 5.0*std::cos(pi*t);
}

int main(int argc, char *argv[]) {

	int params = 2;

	if (argc < params) 
	{
		std::cout << "Por favor, ingrese los parámetros necesarios..." << std::endl;
	} 
	else 
	{
		double xl = atof(argv[1]);
		double xu = atof(argv[2]);

		biseccion<double> bis(xl, xu);
		interpolacion<double> interpol(xl,xu);

		std::cout << "Raíz por bisección en [" << xl << ", " << xu << "]= " << bis(check) << std::endl;
		std::cout << std::endl;
		std::cout << "Raíz por interpolación en [" << xl << ", " << xu << "]= " << interpol(check) << std::endl;
	}

	return 0;
}