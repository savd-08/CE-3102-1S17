#ifndef biseccion_f
#define biseccion_f

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>   
#include <limits>
#include <stdio.h>

/**
*Clase que realiza el metodo de biseccion
*/
template <typename T>
class biseccion { 

private:
	/* Variables utilizadas por el método de bisección */

	T _xl, _xu, _precision;
	int _iters;

	/**
	*Inicializa el funtor del metodo de biseccion
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
	biseccion(T xl, T xu) {

		inicializar(xl, xu);
		
	}

	//Destructor
	~biseccion() {
		
	}

	//Operador del funtor 
	inline T operator() (T (*fnctn)(const T)) {

		std::cout << "i          ";
		std::cout << "xl         ";
		std::cout << "xu         ";
		std::cout << "xr         ";
		std::cout << "ea         " << std::endl;
		std::cout << "_________________________________________________________________" << std::endl;

		T result = calc(fnctn);
		return result;
		
	}


	//funcion para calcular la raiz por metodo de interpolacion lineal
	inline T calc(T (*f)(const T)) {
		
		T raiz = _xl; ///Se inicia con param valido
		T fl = f(_xl); ///sombra 
		T errorAprox = T(); /// error aproximado

		for(int i=_iters ; i>0 ; i--) {
			T raiz_ant(raiz); ///variable para calculo de error
			raiz = (_xl + _xu) / T(2); /// nueva estimacion de la raiz centrada
			T fr = f(raiz); /// sombra de f en el centro

			std::cout << i << "          ";
			std::cout << _xl << "          ";
			std::cout << _xu << "          ";
			std::cout << raiz << "          ";
			std::cout << errorAprox << "          ";

			if (std::abs(raiz) > std::numeric_limits<T>::epsilon()) {
				errorAprox = std::abs((raiz-raiz_ant) / raiz) * T(100); ///nuevo error aprox
			}

			T eval = fl*fr; /// es negativo si extremo inferior
							/// y nuevo centro tienen signos diferentes
			if (eval < T(0)) {
				_xu=raiz;  // es negativo: siga con lado izquierdo
			} else if (eval > T(0)) {
				_xl=raiz;  //es positivo: siga con lado derecho
				fl=fr;
			} else {
				errorAprox = T(0);
				raiz = (std::abs(fl) < std::numeric_limits<T>::epsilon())
				 ? _xl : raiz; //fl == 0
			}

			if (errorAprox < _precision) { return raiz; } //si alcanza la precision estimada, retorne la raiz

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