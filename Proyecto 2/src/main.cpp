#include <iostream>
#include <complex>
#include "rootsalgorithms.h"

using namespace std;
using namespace boost::math::tools;

template <typename T, typename function>
inline void find_roots(polynomial<T> &poly, T* roots, const bool &polish, T x0, function complex_root){
	// grado del polinomio
	int poly_deg = poly.degree();
	polynomial<T> poly_aux(poly);

	for (int i = poly_deg - 1; i >= 0; i--) {
		polynomial<T> p_loop(poly_aux);

		x0 = complex_root(x0); //-----------------------------FUNCION--------------------------


		if (abs(imag(x)) <= 2.0*EPS*abs(real(x))){
			x = Complex(real(x),0.0);
		}

		roots[i]=x;

		poly_aux = deflate(poly_aux, x);

	}

	if (polish){
		for (int i = 0; i < poly_deg; i++){
			laguer(poly,roots[i],its); //-----------------------------FUNCION------------------
		}
	}
}

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

	polynomial<complex<double>> pol{{complex<double>(-82.0,0.0), complex<double>(55.0,0.0), complex<double>(-23.0,0.0),
		complex<double>(-98.0,0.0), complex<double>(112.0,0.0)}};
	polynomial<complex<double>> pol2{{complex<double>(96.0,0.0), complex<double>(-32.0,0.0), complex<double>(45.0,0.0),
		complex<double>(15.0,0.0)}};

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
	polynomial<complex<double>> def_laguerre = deflate(pol2, root_laguerre);
	cout << "Raíces restantes: ";
	for (int i=0; i < def_laguerre.degree() + 1; i++)
			cout << def_laguerre[i] << ", ";
	cout << endl;

	return 0;
}
