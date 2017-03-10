#include "anpi.h"
#include <cmath>

#ifndef ref_f
#define ref_f

namespace anpi{
	namespace ref{

		template <typename T>
		inline T poly_evaluator(T x, T* terms, unsigned int term_amt){
			T result = terms[0] * x;

			for(int i = 1; i < term_amt - 1; i++){
				result = (result + terms[i]) * x;
			}

			result += terms[term_amt - 1];

			return result;
		}

	}
}

#endif
