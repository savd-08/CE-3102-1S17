#include <iostream>
#include <complex>
#include <limits>
#include "rootsalgorithms.h"

using namespace std;
using namespace boost::math::tools;

template <typename T, typename function>
inline void find_roots(polynomial<complex<T>> &poly, complex<T>* roots, const bool &polish, complex<T> x0, function complex_root){
	const double epsilon = std::sqrt(std::numeric_limits<T>::epsilon());
	// grado del polinomio
	int poly_deg = poly.degree();
	polynomial<complex<T>> poly_aux(poly);
	// copia donde se inicia a probar el polinomio
	complex<T> t_root = x0;

	try{
		for (int i = poly_deg-1; i >= 0; i--) {

			t_root = complex_root(x0);

			if (abs(imag(t_root)) <= 2.0*epsilon*abs(real(t_root))){
				t_root = complex<T>(real(t_root),0.0);
			}

			roots[i] = t_root;

			poly_aux = deflate(poly_aux, t_root);

			complex_root = function(poly_aux);
		}
		cout << "Raíces sin pulir: ";
		for(int i = 0; i < 4; i++){
			cout << roots[i] << ", ";
		}
		cout << endl << endl;


	}
	catch(const char* msg){
		cerr << msg << " calculo de raices" << endl;
	}

	try{

		if (polish){
			complex_root = function(poly);
			for (int i = 0; i < poly_deg; i++){
				roots[i] = complex_root(roots[i]);
			}
			cout << "Raíces pulidas:   ";
			for(int i = 0; i < 4; i++){
				cout << roots[i] << ", ";
			}
			cout << endl << endl;
		}

	}
	catch(const char* msg){
		cerr << msg << " pulido" << endl;
	}
}

int main(int argc, char *argv[]) {
	cout.precision(std::numeric_limits<double>::digits10 + 1);
	polynomial<complex<double>> pol1{{complex<double>(1.0,0.0),complex<double>(1.0,0.0), complex<double>(1.0,0.0)}};
	polynomial<complex<double>> pol2{{complex<double>(-105.0,0.0), complex<double>(-1.0,0.0), complex<double>(7.0,0.0),complex<double>(2.0,0.0),complex<double>(1.0,0.0)}};
	polynomial<complex<double>> pol3{{complex<double>(216.0,0.0),complex<double>(-6.0,0.0),complex<double>(-11.0,0.0),complex<double>(1.0,0.0)}};
	/*
	//MULLER
	try{
		muller<double> muller(pol1);
		std::complex<double> root_muller = muller(std::complex<double>(0.0));
		cout << "\n\nRaíz por método de Muller: " << root_muller << "\n";

		laguerre<double> lag(pol1);
		complex<double> root_laguerre = lag(std::complex<double>(0.0));
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
	}   */

	/********************************************************************************/
		//PRUEBA 1

	//MULLER
	//prueba polinomio especificacion
/*	cout << "especificacion" << endl;
	complex<double> *roots1 = new complex<double>[5];
	muller<double>  muller1(pol2);
	find_roots(pol2, roots1, true, complex<double>(0.0), muller1);
	for(int i = 0; i < 5; i++){
		cout << "raices muller: " << roots1[i] << " ";
	}
	cout << endl;

	//LAGUERRE
	cout << "especificacion" << endl;
	complex<double> *roots2 = new complex<double>[5];
	laguerre<double>  laguerre1(pol2);
	find_roots(pol2, roots2, true, complex<double>(0.0), laguerre1);
	for(int i = 0; i < 5; i++){
		cout << "raices laguerre: " << roots2[i] << " ";
	}
	cout << endl;*/

	

	    //PRUEBA 2

	//MULLER
	//prueba polinomio raices reales y complejas
	cout << "raices reales y complejas" << endl;
	complex<double> *roots3 = new complex<double>[4];
	muller<double> muller2(pol2);
	find_roots(pol2, roots3, true, complex<double>(0.0), muller2);
	/*for(int i = 0; i < 4; i++){
		cout << "raices muller: " << roots3[i] << " ";
	}
	cout << endl;

	//LAGUERRE
	//prueba polinomio raices reales y complejas
	cout << "raices reales y complejas" << endl;
	complex<double> *roots4 = new complex<double>[4];
	laguerre<double> laguerre2(pol2);
	find_roots(pol2, roots4, true, complex<double>(0.0), laguerre2);
	for(int i = 0; i < 4; i++){
		cout << "raices laguerre: " << roots4[i] << " ";
	}
	cout << endl;

	

	   //PRUEBA 3

	//MULLER
	//prueba con 3 raíces reales
	cout << "3 raices reales" << endl;
	complex<double> *roots5 = new complex<double>[3];
	muller<double> muller3(pol3);
	find_roots(pol3, roots5, true, complex<double>(0.0), muller3);
	for(int i = 0; i < 3; i++){
		cout << "raices muller: " << roots5[i] << " ";
	}
	cout << endl;   

	//LAGUERRE
	//prueba polinomio raices reales y complejas
	cout << "raices reales y complejas" << endl;
	complex<double> *roots6 = new complex<double>[4];
	laguerre<double> laguerre3(pol3);
	find_roots(pol3, roots6, true, complex<double>(0.0), laguerre3);
	for(int i = 0; i < 4; i++){
		cout << "raices laguerre: " << roots6[i] << " ";
	}
	cout << endl;     */

	return 0;
}
