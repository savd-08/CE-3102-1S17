#include <iostream>
#include <complex>
#include "rootsalgorithms.h"

using namespace std;
using namespace boost::math::tools;

template <typename T, typename function>
inline void find_roots(polynomial<complex<T>> &poly, complex<T>* roots, const bool &polish, complex<T> x0, function complex_root){
	const double epsilon = std::sqrt(std::numeric_limits<T>::epsilon());
	// grado del polinomio
	int poly_deg = poly.degree();
	polynomial<complex<T>> poly_aux(poly);
	// compia donde se inicia a probar el polinomio
	complex<T> t_root = x0;

	try{
		for (int i = poly_deg - 1; i >= 0; i--) {
			polynomial<complex<T>> p_loop(poly_aux);

			t_root = complex_root(x0);

			if (abs(imag(t_root)) <= 2.0*epsilon*abs(real(t_root))){
				t_root = complex<T>(real(t_root),0.0);
			}

			roots[i]=t_root;

			poly_aux = deflate(poly_aux, t_root);

			complex_root = function(poly_aux);

		}
	}
	catch(const char* msg){
		cerr << msg << " calulo de raices" << endl;
	}

	try{

		if (polish){
			for (int i = 0; i < poly_deg; i++){
				complex_root = function(poly);
				complex_root(roots[i]);
			}
		}

	}
	catch(const char* msg){
		cerr << msg << " pulimiento" << endl;
	}
}

int main(int argc, char *argv[]) {
	/*polynomial<double> const a{{-6.0, 1.0, 1.0}};
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
	}*/

	polynomial<complex<double>> pol{{complex<double>(-3.0,0.0), complex<double>(0.0,0.0), complex<double>(1.0,0.0)}};
	polynomial<complex<double>> pol2{{complex<double>(-5.0,0.0), complex<double>(-1.0,0.0), complex<double>(0.0,0.0),complex<double>(2.0,0.0)}};
	polynomial<complex<double>> pol3{{complex<double>(7.0,0.0),complex<double>(0.0,0.0),complex<double>(4.0,0.0),complex<double>(-3.0,0.0),complex<double>(1.0,0.0)}};

	/****************************************************************************************************************/

	//MULLER
	try{
		muller<double> muller(pol);
		std::complex<double> root_muller = muller(std::complex<double>(8.0));
		cout << "\n\nRaíz por método de Muller: " << root_muller << "\n";

		laguerre<double> lag(pol);
		complex<double> root_laguerre = lag(complex<double>(8.0));
		cout << "\n\nRaíz por método de Laguerre: " << root_laguerre << "\n";
	}
	catch(const char* msg){
		cerr << msg << endl;
	}


	try{
		muller<double> muller(pol2);
		std::complex<double> root_muller = muller(std::complex<double>(50.0));
		cout << "\n\nRaíz por método de Muller: " << root_muller << "\n";

		laguerre<double> lag(pol2);
		complex<double> root_laguerre = lag(complex<double>(8.0));
		cout << "\n\nRaíz por método de Laguerre: " << root_laguerre << "\n";
	}
	catch(const char* msg){
		cerr << msg << endl;
	}

	complex<double> *roots = new complex<double>[2];

	muller<double>  muller(pol);

	find_roots(pol, roots, true, complex<double>(10.0), muller);

	for(int i = 0; i < 2; i++){
		cout << "raiz: " << roots[i] << " ";
	}

	cout << endl;

	return 0;
}
