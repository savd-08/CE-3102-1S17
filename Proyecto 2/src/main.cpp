#include "boost_poly.h"
#include <iostream>
#include <complex>
#include "muller.h"

using namespace boost::math::tools;

std::complex<double> f_x(std::complex<double> x){
	return x*x + 3.0;
}

int main(){
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

	anpi::muller<double> muller(f_x);

	std::complex<double> complex1 = muller(0.0, 4.0);

	std::cout <<"\n\nraiz: " << complex1 << "\n";

	return 0;
}
