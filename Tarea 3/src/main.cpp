#include <iostream>
#include "biseccion.h"
#include "interpolacion.h"
#include "newton_raphson.h"


/**************************************************/

//Cálculo de la n-ésima potencia entera de un número de punto flotante.
template <typename fp>
inline fp pow(fp base, int power) {
		fp result = fp(1);
		while(power > 0) {

				if(power & 1) {
						result = (result*base);
				}

				base = (base * base);
				power >>= 1;
		}
		return result;
}

double test_func_1(const double t) {
	static const double pi = 3.14159265358979323856264338327950288;
	return 0.5*std::exp(-t) - 5.0*std::cos(pi*t);
}

double d_test_func_1(const double t) {
	static const double pi = 3.14159265358979323856264338327950288;
	return -0.5*std::exp(-t) + 5.0*pi*std::sin(pi*t);
}

double test_func_2(const double x) {
	return x*x - 3.0;
}

double d_test_func_2(const double x) {
	return 2*x;
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

		biseccion<double> bis(test_func_2);
		interpolacion<double> interpol(test_func_2);
		newton_raphson<double> nt_rp(test_func_2, d_test_func_2);

		std::cout << "Raíz por bisección en [" << xl << ", " << xu << "]= " << bis(xl, xu) << std::endl;
		std::cout << std::endl;
		std::cout << "Raíz por interpolación en [" << xl << ", " << xu << "]= " << interpol(xl, xu) << std::endl;
		std::cout << std::endl;
		std::cout << "Raíz por Newton-Raphson iniciando en " << xl << " = " << nt_rp(xl) << std::endl;
	}

	return 0;
}
