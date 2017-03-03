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
	/**
	 * Evaluación de polinomios por el método de Estrin, utilizando el punto donde
	 * se desea evaluar, los términos del polinomio y la cantidad de dichos términos. 
	 * Se aprovecha la extensión SIMD del procesador. 
	 *  
	 */
	inline double poly_evaluator(double x, double* terms, unsigned int term_amt){
		
		int n_pow2 = nearest_power_2(term_amt);
		int j_limit = n_pow2;
		int num_steps = i_log2(n_pow2);
		double* terms_z = (double*) calloc(n_pow2, sizeof(double));
		__m256d vec_c;
		__m256d vec_x_power = _mm256_set_pd(1, x, 1, x);

		int size_diff = n_pow2 - term_amt; 
		for(int i = 0; i < term_amt; i++){
			terms_z[i+size_diff] = terms[i];
		}
		

		for(int i = 0; i < num_steps; i++){
			j_limit = (j_limit <= 1 ) ? 1 : n_pow2/(1<<(2+i)) ;
			for(int j = 0; j < j_limit; j++){

				vec_c = _mm256_set_pd(terms_z[3 + j*4], terms_z[2 + j*4], terms_z[1 + j*4], terms_z[0 + j*4]);

				vec_c = _mm256_mul_pd(vec_c, vec_x_power);

				vec_c = _mm256_hadd_pd(vec_c, vec_c);

				double* vec_conv = ((double*) (&vec_c));

				terms_z[j*2] = vec_conv[0];
				terms_z[j*2 + 1] = vec_conv[2];
			}
			vec_x_power = _mm256_mul_pd(vec_x_power, vec_x_power);
		}

		x = terms_z[0];
		free(terms_z);

		return x;
	}
	}

	namespace ref{

	inline double poly_evaluator(double x, double* terms, unsigned int term_amt){
		return x;
	}
	}

}
#endif