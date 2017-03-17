#ifndef interpolacion_f
#define interpolacion_f

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>   
#include <limits>
#include <stdio.h>

/**
*Clase que realiza el metodo de interpolacion
*/
template <typename T>
class interpolacion { 

private:
	/* Variables utilizadas por el método de interpolacion */

	T _xl, _xu, _precision;
	int _iters;

	/**
	*Inicializa el funtor del metodo de interpolacion
	*/
	void inicializar(T xl, T xu) {

		/* inicializar variables */
		_xl = xl;
		_xu = xu;
		_precision = sqrt(std::numeric_limits<T>::epsilon());
		_iters = std::numeric_limits<T>::digits;
	}


/*************************************************************/

public:
	//Constructor
	interpolacion(T xl, T xu) {

		inicializar(xl, xu);
		
	}

	//Destructor
	~interpolacion() {
		
	}

	//Operador del funtor 
	inline T operator() (T (*fnctn)(const T)) {

		std::cout << "i          ";
		std::cout << "xl         ";
		std::cout << "xu         ";
		std::cout << "xr         ";
		std::cout << "ea         " << std::endl;
		std::cout << "_________________________________________________________" << std::endl;

		T result = calc(fnctn);
		return result;
	}

		

	//funcion para calcular la raiz por metodo de interpolacion
	inline T calc(T (*f)(const T)) {

		T raiz = _xl; ///Inicio del calculo de la raiz
		T fl = f(_xl); ///funcion evaluada en xl
		T fu = f(_xu); ///funcion evaluada en xu
		T errorAprox; /// error aproximado

		int iu=0, il=0; ///contadores para detectar estancamientos

		for(int i=_iters ; i>0 ; i--) {

			T raiz_ant(raiz); ///variable para calcular el error
			raiz = (_xu-fu*(_xl-_xu))/(fl-fu); /// nueva raiz por interpolacion 
			T fr = f(raiz); /// funcion evaluada en la raiz

			std::cout << i << "          ";
			std::cout << _xl << "          ";
			std::cout << _xu << "          ";
			std::cout << raiz << "          ";
			std::cout << errorAprox << "          ";

			if (std::abs(raiz) > std::numeric_limits<T>::epsilon()) {
				
				errorAprox = std::abs((raiz-raiz_ant)/raiz)*T(100); ///nuevo error aproximado					
			}

			T eval = fl*raiz; 
			if (eval < T(0)) {
				_xu = raiz;
				fu = fr;
				iu = 0;
				il++;
				if (il >= 2) {fl /= T(2);}

			} else if (eval > T(0)) {
				_xl = raiz;
				fl = fr;
				il = 0;
				iu++;
				if(iu >= 2) {fu /= 2;}
			} else {
				errorAprox = T(0);
				raiz = (fl == T(0)) ? _xl : _xu;
			}

			if (errorAprox < _precision) { return raiz; } //si alcanza la precision requerida, retorna raiz

			std::cout << std::endl;
			
		}

		std::cout << 0 << "          ";
		std::cout << _xl << "          ";
		std::cout << _xu << "          ";
		std::cout << raiz << "          ";
		std::cout << errorAprox << "          ";

		return std::numeric_limits<T>::quiet_NaN();
	}

     

};



#endif