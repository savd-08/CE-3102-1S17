#ifndef biseccion_f
#define biseccion_f

#include <cmath>
#include <limits>

/**
*Clase que realiza el metodo de biseccion
*/
template <typename T>
class biseccion { 

private:
	/* Variables utilizadas por el método de bisección */

	T _precision;
	int _iters;
	T (*_func) (const T);


/*************************************************************/

public:
	//Constructor
	biseccion(T (*func) (const T), const T precision = std::sqrt(std::numeric_limits<T>::epsilon()), const int iters = std::numeric_limits<T>::digits) {

		_func = func;
		_precision = precision;
		_iters = iters;
		
	}

	//Destructor
	~biseccion() {
		
	}

	//Operador del funtor 
	inline T operator() (T xl, T xu) {

		T raiz = xl; ///Se inicia con param valido
		T fl = _func(xl); ///sombra 
		T errorAprox = T(); /// error aproximado

		for(int i=_iters ; i>0 ; i--) {
			T raiz_ant(raiz); ///variable para calculo de error
			raiz = (xl + xu) / T(2); /// nueva estimacion de la raiz centrada
			T fr = _func(raiz); /// sombra de f en el centro

			if (std::abs(raiz) > std::numeric_limits<T>::epsilon()) {
				errorAprox = std::abs((raiz-raiz_ant) / raiz) * T(100); ///nuevo error aprox
			}

			T eval = fl*fr; /// es negativo si extremo inferior
							/// y nuevo centro tienen signos diferentes
			if (eval < T(0)) {
				xu=raiz;  // es negativo: siga con lado izquierdo
			} else if (eval > T(0)) {
				xl=raiz;  //es positivo: siga con lado derecho
				fl=fr;
			} else {
				errorAprox = T(0);
				raiz = (std::abs(fl) < std::numeric_limits<T>::epsilon())
				 ? xl : raiz; //fl == 0
			}

			if (errorAprox < _precision) { return raiz; } //si alcanza la precision estimada, retorne la raiz
			
		}

		return std::numeric_limits<T>::quiet_NaN();
	}

};



#endif