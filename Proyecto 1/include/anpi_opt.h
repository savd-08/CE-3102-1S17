#include "anpi.h"
#include "immintrin.h"
#include <cmath>

#ifndef opt_f
#define opt_f

namespace anpi{
	namespace  opt{
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

	}
}

#endif
