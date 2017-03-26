#include <iostream>
#include <complex>
//#include "boost_poly.h"
#include "rootsalgorithms.h"

using namespace std;
using namespace boost::math::tools;

int main(int argc, char *argv[]) {
	polynomial<double> const a{{-6.0, 1.0, 1.0}};
	polynomial<double> const b{{3.0, 1.0, 0.0}};

	polynomial<double> c;
	polynomial<double> d = divide(a, b, c);

	polynomial<double> def = deflate(a, -3.0);


	if(!c.is_zero()){
		std::cout << "Residuo: ";
		for (int i=0; i < c.degree() + 1; i++)
				std::cout << c[i] << ", ";
		std::cout << std::endl;
	}

	if(!d.is_zero()){
		std::cout << "Cociente div: ";
		for (int i=0; i < d.degree() + 1; i++)
				std::cout << d[i] << ", ";
		std::cout << std::endl;
	}

	if(!def.is_zero()){
		std::cout << "Cociente def: ";
		for (int i=0; i < def.degree() + 1; i++)
				std::cout << def[i] << ", ";
		std::cout << std::endl;
	}

	polynomial<std::complex<double>> pol{{std::complex<double>(-3.0,0.0), std::complex<double>(0.0,0.0), std::complex<double>(1.0,0.0)}};
	polynomial<std::complex<double>> pol2{{std::complex<double>(6.0,0.0), std::complex<double>(4.0,0.0), std::complex<double>(6.0,0.0), std::complex<double>(8.0,0.0)}};

	/****************************************************************************************************************/

	int params = 1;
	if (argc < params) {cout << "Por favor, ingrese los parámetros necesarios..." << endl;}
	else
	{
		//MULLER
		muller<double> muller(pol2);
		polynomial<complex<double>> root_mul = muller(0.0);
		std::cout << "\n\nRaíz por método de Muller: " << root_mul[0] << " " << root_mul[1] << "\n";

		//LAGUERRE
		bool complex_mode = (char)atoi(argv[1]);
		laguerre<double> lag(pol2, complex_mode);
		polynomial<complex<double>> root_lag = lag(0.0);
		cout << "\n\nRaíz por método de Laguerre: " << root_lag[0] << " " << root_lag[1] << "\n";

	}

	return 0;
}
