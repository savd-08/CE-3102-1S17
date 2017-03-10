#include <cmath>

#ifndef anpi_f
#define anpi_f

namespace anpi{

	//Encuentra la potencia de 2 más cercana hacia arriba de un entero.
	//Obtenido de http://graphics.stanford.edu/~seander/bithacks.html
	inline unsigned int nearest_power_2 (unsigned int val){

		val--;
		val |= val >> 1;
		val |= val >> 2;
		val |= val >> 4;
		val |= val >> 8;
		val |= val >> 16;
		val++;

		return val;
	}


	//Calculo del logaritmo base 2 de un número entero.
	//Obtenido de http://graphics.stanford.edu/~seander/bithacks.html
	inline unsigned int i_log2(unsigned int val){
		unsigned int r;
		static const int MultiplyDeBruijnBitPosition2[32] =
		{
		  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
		  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
		};
		r = MultiplyDeBruijnBitPosition2[(unsigned int)(val * 0x077CB531U) >> 27];
		return r;
	}


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

	//Cálculo del factorial
	template<typename T>
    T factorial(const unsigned int n) {
        return n < 2 ? T(1) : static_cast<T>(n) * factorial<T>(n-1);
    }

	/**
	* Functor que calcula los terminos de taylor para el logartimo.
	* Llama el algoritmo de Horner para el calculo del logaritmo.
	*/
	template <class T>
	class ln_a {
		private:
		//arreglo con los monomios
		T* _coef;
		//cantidad de terminos de la serie
		unsigned int _terms;
		//centro de la serie
		T _center;
		//puntero a la funcion evaluadora
		T (*_poly_evaluator)(T, T*, unsigned int);

		//metodo inicializador
		void init(const T center, unsigned int terms, T(*poly_evaluator)(T, T*, unsigned int)){
			_center = center;
			_terms = terms;
			_coef = new T[_terms];
			_poly_evaluator = poly_evaluator;

			//calulo de los coeficientes
			for(int i = (_terms - 1); i >= 0; i--){
				_coef[(_terms - 1) - i] = diff(_center, i);
			}
		}

		public:
			//constructor
			ln_a(T center, int terms, T(*poly_evaluator)(T, T*, unsigned int)){
				init(center, terms, poly_evaluator);
			}

			/// Destructor
					~ln_a() { free(_coef); }

			//functor
			inline T operator()(T value){
				//valor que va elevado
				T h = value - _center;
				return (*_poly_evaluator)(h, _coef, _terms);
			}//termina sobrecarga

			//enesima derivada
			inline T diff(T x, int n){
				T result = 0;
				//Para el caso de n == 0
				if(n == 0){
					result = std::log(x);
				}
				//de 0 en adelante
				else{
					result = 1 / (n * anpi::pow(x,n));
					result = (!(n & 1)) ? -result : result;
				}
				return result;
				}//termina funcion diff
	};//fin del functor del logaritmo

	/**
	* Functor que calcula los terminos de taylor para el coseno.
	* Llama el algoritmo de Horner para el calculo del coseno.
	*/
	template <typename T>
	class cos_a {

		private:
			/// Coeficientes de la serie de Taylor
			T* _coef;

			/// Centro de la serie
			T _center;
			/// Cantidad de terminos de la serie
			unsigned int _terms;

			/// Puntero a la la funcion evaluadora
			T (*_poly_evaluator)(T, T*, unsigned int);

			/**
				*Inicializa los coeficientes de la serie de Taylor centrada en center
				* para el numero dado de términos
				*/
			void init(const T center, unsigned int terms, T(*poly_evaluator)(T, T*, unsigned int)) {

					_center = center;
					_terms = terms;
					_coef = new T[_terms];
					_poly_evaluator = poly_evaluator;
					//Solicita memoria para la creacion de los coeficientes


					//Calculo de los coeficientes
					for (int i = (terms - 1); i >= 0; i--) {
							_coef[(terms - 1) - i] = (diff(center, i)) / (factorial<T>(i));

					}

			}

		public:
			/**
			* Unico constructor obliga a dar centro y numero de terminos
			*  de la aproximacion
			*/
			cos_a(const T center, int terms, T(*poly_evaluator)(T, T*, unsigned int)) {
					//Se inicializa la clase
					init(center, terms, poly_evaluator);
			}

			/// Destructor
			~cos_a() { free(_coef); }

			/// Evaluación de la función coseno
			inline T operator() (T val) {
					T h = val - _center;
					return (*_poly_evaluator)(h, _coef, _terms);
			}

			/// Evaluación de la n-ésima derivada de coseno
			inline T diff(T x, int n) {

					T result;
					if (n & 1) {
						unsigned int parity = (n + 1) >> 1;
						result = (parity & 1) ?  -std::sin(x) : std::sin(x);
					} else {

						unsigned int parity = n >> 1;
						result = (parity & 1) ?  -std::cos(x) : std::cos(x);
					}

					return result;
			}
	};

}
#endif
