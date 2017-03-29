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

	polynomial<complex<double>> pol1{{complex<double>(-8.0,0.0), complex<double>(-4.0,0.0), complex<double>(2.0,0.0), complex<double>(-1.0,0.0), complex<double>(1.0,0.0)}};
	polynomial<complex<double>> pol2{{complex<double>(-105.0,0.0), complex<double>(-1.0,0.0), complex<double>(7.0,0.0),complex<double>(2.0,0.0),complex<double>(1.0,0.0)}};
	polynomial<complex<double>> pol3{{complex<double>(216.0,0.0),complex<double>(-6.0,0.0),complex<double>(-11.0,0.0),complex<double>(1.0,0.0)}};

	/****************************************************************************************************************/

	//MULLER
	try{
		muller<double> muller(pol1);
		std::complex<double> root_muller = muller(std::complex<double>(8.0));
		cout << "\n\nRaíz por método de Muller: " << root_muller << "\n";

		laguerre<double> lag(pol1);
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


	//prueba polinomio especificacion
	cout << "especificaion" << endl;
	complex<double> *roots1 = new complex<double>[4];
	muller<double>  muller1(pol1);
	find_roots(pol1, roots1, true, complex<double>(0.0), muller1);
	for(int i = 0; i < 4; i++){
		cout << "raiz: " << roots1[i] << " ";
	}
	cout << endl;

	//prueba polinomio raices reales y complejas
	cout << "raices reales y complejas" << endl;
	complex<double> *roots2 = new complex<double>[4];
	muller<double> muller2(pol2);
	find_roots(pol2, roots2, true, complex<double>(0.0), muller2);
	for(int i = 0; i < 4; i++){
		cout << "raiz: " << roots2[i] << " ";
	}
	cout << endl;

	//prueba con 3 raíces reales
	cout << "3 raices reales" << endl;
	complex<double> *roots3 = new complex<double>[3];
	muller<double> muller3(pol3);
	find_roots(pol3, roots3, true, complex<double>(0.0), muller3);
	for(int i = 0; i < 3; i++){
		cout << "raiz: " << roots3[i] << " ";
	}
	cout << endl;

	return 0;
}
