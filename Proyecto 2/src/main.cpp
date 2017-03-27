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
	muller<double> muller(pol);
	std::complex<double> root_muller = muller(0.0);
	cout << "\n\nRaíz por método de Muller: " << root_muller << "\n";

	//Deflacion para calcular las demas raices 
	polynomial<complex<double>> def_muller = deflate(pol, root_muller);
	cout << "Raíces restantes: ";
	for (int i=0; i < def_muller.degree() + 1; i++)
			cout << def_muller[i] << ", ";
	cout << endl;

	/*********************************************************/  

	//LAGUERRE
	laguerre<double> lag(pol);
	complex<double> root_laguerre = lag(0.0);
	cout << "\n\nRaíz por método de Laguerre: " << lag(0.0) << "\n";

	//Deflacion para calcular las demas raices 
	polynomial<complex<double>> def_laguerre = deflate(pol, root_laguerre);
	cout << "Raíces restantes: ";
	for (int i=0; i < def_laguerre.degree() + 1; i++)
			cout << def_laguerre[i] << ", ";
	cout << endl;



	return 0;
}
