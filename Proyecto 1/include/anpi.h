#include "immintrin.h"

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

	namespace opt{
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
				//std::cout << "n_pow2: "<< n_pow2 << std::endl;
				//int n_pow2_cp = n_pow2;

				int zeros = n_pow2 - term_amt;
				//std::cout << "Zeros: "<< zeros << std::endl;

				int num_steps = i_log2(n_pow2);
				//std::cout << "num_steps: "<< num_steps << std::endl;

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
					//std::cout << "Ignore: "<< ignore << std::endl;
					int j = 0;
					do{
						//std::cout << "........................j=: "<< j << std::endl;
						int j_step = (zeros+j)<<2;
						coef_vector.set(terms_z[3 + j_step], terms_z[2 + j_step], terms_z[1 + j_step], terms_z[j_step]);

						/*std::cout << "coef_vector: ";
							for(int k = 0; k < 4; k++){
								std::cout << coef_vector.data[k] << ", ";
							}
						std::cout << std::endl;*/

						coef_vector.mul(coef_vector.data, x_power_vector.data);

						coef_vector.hadd(coef_vector.data, coef_vector.data);

						terms_z[(zeros<<1) + (j<<1)] = coef_vector.data[0];
						terms_z[(zeros<<1) + (j<<1) + 1] = coef_vector.data[type_index];

						/*std::cout << "terms_z: ";
							for(int k = 0; k < n_pow2; k++){
								std::cout << terms_z[k] << ", ";
							}
						std::cout << std::endl;*/

						j++;
					} while(j < ((n_pow2)>>(2+i)) - zeros);
					//n_pow2_cp >>= 1; 
					zeros >>= 1;
					x_power_vector.mul(x_power_vector.data, x_power_vector.data);
				}

				x = terms_z[0];
				free(terms_z);

				return x;
		}
	}

	namespace ref{

		template <typename T>
		inline T poly_evaluator(T x, T* terms, unsigned int term_amt){
			T result = terms[0]*x;
			for(int i = 1; i < term_amt-1; i++){
				result = (result + terms[i])*x;
			}
			result += terms[term_amt-1];
			return result;
		}
	}

}
#endif