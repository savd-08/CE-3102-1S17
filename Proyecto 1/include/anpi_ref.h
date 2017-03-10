#include "anpi.h"
#include <cmath>

#ifndef ref_f
#define ref_f

namespace anpi::ref{

	template <typename T>
	inline T poly_evaluator(T x, T* terms, unsigned int term_amt){
		T result = terms[0] * x;

		for(int i = 1; i < term_amt - 1; i++){
			result = (result + terms[i]) * x;
		}

		result += terms[term_amt - 1];

		return result;
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

		//metodo inicializador
		void init(const T center, unsigned int terms){
			_center = center;

			_terms = terms;

			_coef = new T[_terms];

			//calulo de los coeficientes
			for(int i = (_terms - 1); i >= 0; i--){
				_coef[(_terms - 1) - i] = diff(_center, i);
			}
		}

	public:
		//constructor
		ln_a(T center, int terms){
			init(center, terms);
		}

		/// Destructor
        ~ln_a() { free(_coef); }

		//functor
		inline T operator()(T value){
			//valor que va elevado
			T h = value - _center;
			return anpi::ref::poly_evaluator(h, _coef, _terms);

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

        /**
          *Inicializa los coeficientes de la serie de Taylor centrada en center
          * para el numero dado de términos
          */
        void init(const T center, unsigned int terms) {

            _center = center;
            _terms = terms;
            _coef = new T[_terms];

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
        cos_a(const T center, int terms) {
            //Se inicializa la clase
            init(center, terms);
        }

        /// Destructor
        ~cos_a() { free(_coef); }

        /// Evaluación de la función coseno
        inline T operator() (T val) {
            T h = val - _center;
            return anpi::ref::poly_evaluator(h, _coef, _terms);
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