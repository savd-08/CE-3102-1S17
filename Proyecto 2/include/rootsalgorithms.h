#include <iostream>
#include <complex>
#include <limits>
#include <cmath>
#include <fstream>
#include "boost_poly.h"


namespace boost{ namespace math{ namespace tools{

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

	//Genera archivo con los puntos para las gráficas
	template <typename T>
	void load_plot(T y, T x, std::ofstream& file) {

		//Introducir datos a los archivos
		file << (x) << "\t";
		file << (y) << std::endl;
	}

/********************************************************************************************************************************/

	/**
	*Clase que realiza el metodo de Laguerre
	*/
	template <class T>
	class laguerre {

	private:
		polynomial<std::complex<T>> _poly;
		T _precision; //Tolerancia
		unsigned int _iters; //Iteraciones
		std::complex<T>* d_coefs; //Arreglo para almacenar coeficientes de la primera derivada

	public:
		//Constructor
		laguerre(polynomial<std::complex<T>> ppoly, unsigned int piters = 100000, T pprecision = std::sqrt(std::numeric_limits<T>::epsilon())) {
			//Se establecen parametros
			_iters = piters;
			_precision = pprecision;
			_poly = ppoly;
			d_coefs = new std::complex<T>[_poly.degree()];
		}

		//Destructor
		~laguerre() {
			//free(d_coefs);
		}

		/***************FUNCIONES PARA EVALUAR POLINOMIOS***************************/

			//Evalua la funcion derivada en un punto especifico
		std::complex<T> dpol(std::complex<T> x, const polynomial<std::complex<T>>& ppol) {

			std::complex<T> result = std::complex<T>(0);
			for (int i = ppol.degree(); i > 0; i--) {
				result += ppol[i] * std::complex<T>(i) * pow(x, i-1); //Calculo de la primera derivada del polinomio
				d_coefs[i-1] = ppol[i] * std::complex<T>(i); //Almacena los coeficientes para utilizarlos en 'd2pol(T x)'
			}
			return result;
		}

		//Evalua la funcion derivada dos veces en un punto especifico
		std::complex<T> d2pol(std::complex<T> x) {

			int _length = (sizeof(d_coefs)/sizeof(*d_coefs));
			std::complex<T> result = std::complex<T>(0);
			for (int i = _length+1; i > 0; i--) {
				result += d_coefs[i] * std::complex<T>(i) * pow(x, i-1); //Calculo de la segunda derivada del polinomio
			}
			return result;
		}
		/***************************CALCULO DE LA RAIZ**********************************/

		//Operador del funtor
		std::complex<T> operator() (std::complex<T> lim_inf) {

			/*****Método de Laguerre para la busqueda de raices*****/
			std::complex<T> tmp_root, new_root, g, h, c, error; //Variables para calcular la aproximacion de la raiz
			tmp_root = lim_inf;   //Se inicia con una raiz estimada igual al limite inferior del rango
			std::complex<T> _degree = std::complex<T>(_poly.degree()); //Grado del polinomio a evaluar
			std::ofstream tmp_file("laguerre.txt"); //Archivo para la graficacion

			//Iteraciones para buscar la raiz mas aproximada dentro del determinado rango
			for (int i = 0; i < _iters; i++) {

				if(i >= 1) g = std::complex<T>(-1) * dpol(tmp_root, _poly) / poly_evaluator(tmp_root, _poly); //Verificacion de signo para
				else g = dpol(tmp_root, _poly) / poly_evaluator(tmp_root, _poly);              // la variable G

				h = pow(g,2) - (d2pol(tmp_root) / poly_evaluator(tmp_root, _poly)); //calculo del termino H para el metodo Laguerre

				/*Se evalua un termino que compone la variable c para determinar
				el numerador mas grande en la operación */
				std::complex<T> term1 = g + sqrt((_degree-std::complex<T>(1))*((_degree*h)-(pow(g,2))));
				std::complex<T> term2 = g - sqrt((_degree-std::complex<T>(1))*((_degree*h)-(pow(g,2))));
				if (std::abs(term1) > std::abs(term2)) c = term1 / _degree;
				else c = term2 / _degree;

				new_root = tmp_root + (std::complex<T>(1) / c); //Se establece la nueva raiz aproximada
				error = std::abs((new_root - tmp_root) / new_root);
				error = poly_evaluator(new_root, _poly);

				load_plot(new_root.imag(), new_root.real(), tmp_file); //Carga los datos calculados en un archivos para las graficas

				// Se verifica si ya la raiz converge lo suficiente con un valor estimado establecido
				if ( std::abs(error) <= _precision ) {
					if(std::abs(std::real(new_root)) <= _precision){
		      			new_root = std::complex<T>(0.0, std::imag(new_root));
		      		}
		      		if(std::abs(std::imag(new_root)) <= _precision){
		      			new_root = std::complex<T>(std::real(new_root), 0.0);
		    		}

		      		tmp_file.close(); //Cierra el archivo
					return new_root; //Raiz aproximada por metodo de Laguerre

				 //la raiz temporal se cambia por la nueva raiz, para repetir el proceso
				}

				tmp_root = new_root;
			}

			throw("Cantidad maxima de iteraciones alcanzada");
   		}//Fin del operador
	}; //laguerre

/********************************************************************************************************************************/

	  template <class T>
	  class muller
	  {
	    private:
	      //limite de iteraciones
	      unsigned int _iters;
	      //tolerancia
	      T _precision;
	      //polinomio
	      polynomial<std::complex<T>> _polinomio;

	    public:

	      //constructor
	      muller(polynomial<std::complex<T>> polinomio, unsigned int iters = 100000, T precision = std::sqrt(std::numeric_limits<T>::epsilon())){
					_precision = precision;
					_iters = iters;
					_polinomio = polinomio;
	      }

	      std::complex<T> operator()(std::complex<T> x1){

					T diff = T(0.001);

					//paso a numeros complejos
					std::complex<T> cx0 = x1-diff;
					std::complex<T> cx1 = x1;
					std::complex<T> cx2 = x1+diff;

					std::ofstream tmp_file("muller.txt"); //Archivo para la graficacion

					for(unsigned int i = 2; i < _iters; i++){
					  //calculo de distancias
					  std::complex<T> h1 = cx1 - cx0;
					  std::complex<T> h2 = cx2 - cx1;
					  //cambio de la funcion
					  std::complex<T> delta1 = (poly_evaluator(cx1, _polinomio) - poly_evaluator(cx0, _polinomio))/h1;
					  std::complex<T> delta2 = (poly_evaluator(cx2, _polinomio) - poly_evaluator(cx1, _polinomio))/h2;
					  //paso
					  std::complex<T> d = (delta2 - delta1)/(h2+h1);
					  //calculo de b
					  std::complex<T> b = delta2 + (h2*d);
					  //calculo de la raíz cuadrada del denomidador
					  std::complex<T> D = std::sqrt((b*b)-(4.0*poly_evaluator(cx2, _polinomio)*d));
					  //temporales para validación de signo
					  std::complex<T> tmp1 = b - D;
					  std::complex<T> tmp2 = b + d;
					  //denominador de la resta
					  std::complex<T> E = 1;

					  if(std::abs(tmp1) < std::abs(tmp2) ){
					    E = b + D;
					  }
					  else{
					    E = b - D;
					  }
					  //parte que se le resta a x2
					  std::complex<T> h = ((-2.0) * poly_evaluator(cx2, _polinomio))/E;
					  //nueva x
					  std::complex<T> p = cx2 + h;

					  load_plot(p.imag(), p.real(), tmp_file); //Carga los datos calculados en un archivos para las graficas

						std::complex<T> y_p = poly_evaluator(p, _polinomio);

					  if(std::abs(y_p) <= _precision){
					    if(std::abs(std::real(p)) <= _precision){
					      p = std::complex<T>(0.0, std::imag(p));
					    }
					    if(std::abs(std::imag(p)) <= _precision){
					      p = std::complex<T>(std::real(p), 0.0);
					    }

					    tmp_file.close();
					    return p;
					  }

					  cx0 = cx1;
					  cx1 = cx2;
					  cx2 = p;

					}

					throw("Cantidad maxima de iteraciones alcanzada");
	      }
		};//fin de muller

} //namespace tools
} //namespace math
} //namespace boost
