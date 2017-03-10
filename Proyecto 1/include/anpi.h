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

}
#endif
