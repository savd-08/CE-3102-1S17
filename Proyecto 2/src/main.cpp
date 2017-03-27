#include <iostream>
#include <complex>
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
		cout << "Residuo: ";
		for (int i=0; i < c.degree() + 1; i++)
				cout << c[i] << ", ";
		cout << endl;
	}

	if(!d.is_zero()){
		cout << "Cociente div: ";
		for (int i=0; i < d.degree() + 1; i++)
				cout << d[i] << ", ";
		cout << std::endl;
	}

	if(!def.is_zero()){
		cout << "Cociente def: ";
		for (int i=0; i < def.degree() + 1; i++)
				cout << def[i] << ", ";
		cout << std::endl;
	}

	polynomial<complex<double>> pol{{complex<double>(-3.0,0.0), complex<double>(0.0,0.0), complex<double>(1.0,0.0)}};
	polynomial<complex<double>> pol2{{complex<double>(6.0,0.0), complex<double>(4.0,0.0), complex<double>(6.0,0.0), complex<double>(8.0,0.0)}};

	/****************************************************************************************************************/

	//MULLER
	muller<double> muller(pol2);
	cout << "\n\nRaíz por método de Muller: " << muller(0.0) << "\n";

	//LAGUERRE
	laguerre<double> lag(pol2);
	cout << "\n\nRaíz por método de Laguerre: " << lag(0.0) << "\n";


	return 0;
}
