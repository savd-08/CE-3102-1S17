#ifndef interpolacion_f
#define interpolacion_f

#include <cmath>
#include <limits>
#include <iostream>

/**
*Clase que realiza el metodo de interpolacion
*/
template <typename T>
class interpolacion { 

private:
	/* Variables utilizadas por el método de interpolacion */

	T _precision;
	int _iters;
	T (*_func) (const T);

/*************************************************************/

public:
	//Constructor
	interpolacion(T (*func) (const T), const T precision = std::sqrt(std::numeric_limits<T>::epsilon()), const int iters = std::numeric_limits<T>::digits) {

		_func = func;
		_precision = precision;
		_iters = iters;
		
	}

	//Destructor
	~interpolacion() {
		
	}

	//Operador del funtor 
	T operator() (T xl, T xu) {

		T raiz = xl; ///Inicio del calculo de la raiz
		T fl = _func(xl); ///funcion evaluada en xl
		T fu = _func(xu); ///funcion evaluada en xu
		T errorAprox; /// error aproximado

		int iu=0, il=0; ///contadores para detectar estancamientos

		for(int i=_iters ; i>0 ; i--) {

			T raiz_ant(raiz); ///variable para calcular el error
			raiz = xu - fu * (xl - xu) / (fl - fu) ; /// nueva raiz por interpolacion 
			T fr = _func(raiz); /// funcion evaluada en la raiz

			if (std::abs(raiz) > std::numeric_limits<T>::epsilon()) {
				
				errorAprox = std::abs((raiz-raiz_ant)/raiz)*T(100); ///nuevo error aproximado
			}

			T eval = fl*raiz; 
			if (eval < T(0)) {
				xu = raiz;
				fu = fr;
				iu = 0;
				il++;
				if (il >= 2) {fu /= T(2);}
			} else if (eval > T(0)) {
				xl = raiz;
				fl = fr;
				il = 0;
				iu++;
				if(iu >= 2) {fl /= T(2);}
			} else {
				errorAprox = T(0);
				raiz = (std::abs(fl) < std::numeric_limits<T>::epsilon()) ? xl : xu;
			}

			if (errorAprox < _precision) { 
				return raiz; 
			} //si alcanza la precision requerida, retorna raiz
			
		}

		return std::numeric_limits<T>::quiet_NaN();
	}     

};

#endif