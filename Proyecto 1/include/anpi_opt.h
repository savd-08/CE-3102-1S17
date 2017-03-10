#include "anpi.h"
#include "immintrin.h"
#include <cmath>

#ifndef opt_f
#define opt_f

namespace anpi::opt{
	template<typename T>
	class simd;

	template<>
	struct simd<double>
	{
	    typedef __m256d simd_vector;
	    simd_vector data;

	    void set(const double &a, const double &b, const double &c, const double &d){
	       data = _mm256_set_pd(a, b, c, d);
	    }

	    void add(const simd_vector &a, const simd_vector &b){
	       data = _mm256_add_pd(a, b);
	    }

	    void mul(const simd_vector &a, const simd_vector &b){
	       data = _mm256_mul_pd(a, b);
	    }

	    void hadd(const simd_vector &a, const simd_vector &b){
	       data = _mm256_hadd_pd(a, b);
	    }
	};

	template<>
	struct simd<float>
	{
	    typedef __m128 simd_vector;
	    simd_vector data;


	    void set(const float &a, const float &b, const float &c, const float &d){
	       data = _mm_set_ps(a, b, c, d);
	    }

	    void add(const simd_vector &a, const simd_vector &b){
	       data = _mm_add_ps(a, b);
	    }

	    void mul(const simd_vector &a, const simd_vector &b){
	       data = _mm_mul_ps(a, b);
	    }

	    void hadd(const simd_vector &a, const simd_vector &b){
	       data = _mm_hadd_ps(a, b);
	    }
	};


	/**
	 * Evaluación de polinomios por el método de Estrin, utilizando el punto donde
	 * se desea evaluar, los términos del polinomio y la cantidad de dichos términos.
	 * Se aprovecha la extensión SIMD del procesador.
	 *
	 */
	template <typename T>
	inline T  poly_evaluator(T x, T* terms, unsigned int term_amt){

		int n_pow2 = nearest_power_2(term_amt);

		int zeros = n_pow2 - term_amt;

		int num_steps = i_log2(n_pow2);

		T* terms_z = (T*) calloc(n_pow2, sizeof(T));

		constexpr unsigned int type_index = sizeof(T)/sizeof(float);

		simd<T> coef_vector;
		simd<T> x_power_vector;
		x_power_vector.set(T(1), x, T(1), x);

		for(int i = 0; i < term_amt; i++){
			terms_z[i+zeros] = terms[i];
		}

		zeros >>= 2;

		for(int i = 0; i < num_steps; i++){
			int j = 0;
			do{
				int j_step = (zeros+j)<<2;

				coef_vector.set(terms_z[3 + j_step], terms_z[2 + j_step], terms_z[1 + j_step], terms_z[j_step]);

				coef_vector.mul(coef_vector.data, x_power_vector.data);

				coef_vector.hadd(coef_vector.data, coef_vector.data);

				terms_z[(zeros<<1) + (j<<1)] = coef_vector.data[0];
				terms_z[(zeros<<1) + (j<<1) + 1] = coef_vector.data[type_index];

				j++;

			} while(j < ((n_pow2)>>(2+i)) - zeros);

			zeros >>= 1;
			x_power_vector.mul(x_power_vector.data, x_power_vector.data);
		}

		x = terms_z[0];
		free(terms_z);

		return x;
	}

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
			return anpi::opt::poly_evaluator(h, _coef, _terms);

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
    * Llama el algoritmo de Estrin optimizado para el calculo del coseno.
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
            return anpi::opt::poly_evaluator(h, _coef, _terms);
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