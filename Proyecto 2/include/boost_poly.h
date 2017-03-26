#include <boost/math/tools/polynomial.hpp>
#include <iostream>

namespace boost{ namespace math{ namespace tools{
	/*
* Algoritmo general de división polinomial. Recibe el dividendo (n), el divisior (d)
* y el residuo (r) de la divisón como argumentos. Retorna el cociente.
* Todos los argumentos son del tipo boost::math::tools::polynomial.
*/
template <typename T>
inline polynomial<T> divide(const polynomial<T>& n, const polynomial<T>& d, polynomial<T>& r){

	//Se obtienen los grados de cada polinomio
	int n_deg = n.degree();
	int d_deg = d.degree();
	int q_deg = n_deg - d_deg;

	//Se incializan contadores
	int ds_deg;
	int i;

	//Se copia el dividendo
	r = polynomial<T>(n);

	//Instanciación de cociente como copia de n
	polynomial<T> q(n);

	//Instanciación del polinomio divisor auxiliar como copia de n
	polynomial<T> d_shifted(n);

	//Se inicializa el cociente en cero
	for (i = 0; i <= n_deg; i++){
		q[i] = T(0);
	}

	while(n_deg >= d_deg){

		//Se se establece el divisor auxiliar en cero
		for (i = 0; i <= d_deg; i++){
			d_shifted[i] = T(0);
		}

		//Se realiza una copia del divisor, desplazada (n_deg-d_deg) veces
		//Equivalente a multiplicar por x^(n_deg-d_deg)
		for(i = 0; i <= d_deg; i++){
			d_shifted[i + n_deg - d_deg] = d[i];
		}


		ds_deg = n_deg;

		//Cálculo del elemento de mayor orden del cociente
		q[n_deg - d_deg] = r[n_deg]/d_shifted[ds_deg];

		//Multiplicación de todo el divisor auxiliar por el cociente
		for(i = 0; i <= q_deg ; i++) {
			d_shifted[i] = d_shifted[i] * q[n_deg - d_deg];
		}

		//Se resta al residuo parcial el divisor auxiliar, lo que elimina el termino de mayor orden
		for(i = 0 ; i <= n_deg; i++) {
			r[i] = r[i] - d_shifted[i];
		}

		n_deg--;
	}

	//Se eliminan los coeficientes iguales a cero
	r.normalize();
	q.normalize();
	return q;
}


/*
* Función para deflación de polinomios. Se toma la deflación como un caso
* especial de división polinomial con divisor de grado 1 y residuo 0.
* Recibe el dividendo (n) y el divisior de grado 1 o raíz (d) como argumentos.
* Retorna el cociente. Todos los polinomios se representan con boost::math::tools::polynomial.
*/
template <typename T>
inline polynomial<T> deflate(const polynomial<T>& n, const T d){
	//Se obtiene el grado del dividendo
	int n_deg = n.degree();

	//Se genera el cociente como una copia del dividendo
	polynomial<T> q(n);

	//Coeficiente que debe restarse al dividendo en cada iteración
	T r = q[n_deg];

	//Se establece el cociente con un grado menor que el dividendo
	q[n_deg] = T(0);

	T c;

	for(int i = n_deg-1; i >= 0; i--) {
		//Valor del coeficiente de i-ésimo grado
		c = q[i];

		//Coeficiente del cociente
		q[i] = r;

		//Residuo de la división polinomial
		r = d*r + c;
	}

	//Se eliminan los coeficientes iguales a cero
	q.normalize();
	return q;
}

/**
* Algoritmo de horner para evaluar polinomios
*/
template <typename T>
inline T poly_evaluator(T x, const polynomial<T>& terms){

	int term_amt = terms.degree()+1;

	T result = terms[term_amt - 1] * x;

	for(int i = term_amt - 2; i > 0; i--){
		result = (result + terms[i]) * x;
	}

	result += terms[0];

	return result;
}


} //namespace tools
} //namespace math
} //namespace boost
