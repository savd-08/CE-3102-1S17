#ifndef nr_f
#define nr_f

#include <cmath>
#include <limits>

/**
*Clase que realiza el metodo de Newton-Raphs
*/
template <typename T>
class newton_raphson { 

private:
	/* Variables utilizadas por el método de interpolacion */

	T _precision;
	int _iters;
	T (*_func) (const T);
	T (*_d_func) (const T);

/*************************************************************/

public:
	//Constructor
	newton_raphson(T (*func) (const T), T (*d_func) (const T), const T precision = std::sqrt(std::numeric_limits<T>::epsilon()), const int iters = std::numeric_limits<T>::digits) {

		_func = func;
		_d_func = d_func;
		_precision = precision;
		_iters = iters;
		
	}

	//Destructor
	~newton_raphson() {
		
	}

	//Operador del funtor 
	T operator() (T xs) {

		T x0 = (xs == 0) ? std::numeric_limits<T>::epsilon()*1000 : xs; ///Inicio del calculo de la raiz
		T x1 = T(1);
		T errorAprox; /// error aproximado

		for(int i=_iters ; i > 0 ; i--) {
			
			T y = _func(x0);
			T dy = _d_func(x0);

			if (std::abs(dy) < std::numeric_limits<T>::epsilon()) {
				return std::numeric_limits<T>::quiet_NaN();				
			}

			x1 = x0 - (y/dy);

			if (std::abs(x1) > std::numeric_limits<T>::epsilon()) {
				errorAprox = std::abs((x1-x0)/x1)*T(100); ///nuevo error aproximado
			}

			if (errorAprox < _precision) { return x1; } //si alcanza la precision requerida, retorna raiz

			x0 = x1;
			
		}

		return std::numeric_limits<T>::quiet_NaN();
	}     

};

#endif