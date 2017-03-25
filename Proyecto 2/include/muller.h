#include <complex>
#include <limits>

#ifndef anpi_f
#define anpi_f

namespace anpi{

  template <class T>
  class muller
  {
    private:
      //limite de iteraciones
      unsigned int _iters;
      //tolerancia
      T _precision;
      //funcion
      std::complex<T>(*f_x)(std::complex<T> x);

    public:

      //constructor
      muller(std::complex<T>(*pf_x)(std::complex<T> x),T precision = std::sqrt(std::numeric_limits<T>::epsilon()), unsigned int iters = 1000){
        _precision = precision;
        _iters = iters;
        f_x = pf_x;
      }

      std::complex<T> operator()(T x1, T x2){

        //paso a numeros complejos
        std::complex<T> cx0 = x1;
        std::complex<T> cx2 = x2;

        std::complex<T> cy0 = f_x(cx0);
        std::complex<T> cy2 = f_x(cx2);

        //secante para calcular un tercer punto
        std::complex<T> cx1 = (cx2 - cy2)*((cx2 - cx0)/(cy2 - cy0));

        for(unsigned int i = 2; i < _iters; i++){
          //calculo de distancias
          std::complex<T> h1 = cx1 - cx0;
          std::complex<T> h2 = cx2 - cx1;
          //cambio de la funcion
          std::complex<T> delta1 = (f_x(cx1) - f_x(cx0))/h1;
          std::complex<T> delta2 = (f_x(cx2) - f_x(cx1))/h2;
          //paso
          std::complex<T> d = (delta2 - delta1)/(h2+h1);
          //calculo de b
          std::complex<T> b = delta2 + (h2*d);
          //calculo de la raíz cuadrada del denomidador
          std::complex<T> D = std::sqrt((b*b)-(4.0*f_x(cx2)*d));
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
          std::complex<T> h = ((-2.0) * f_x(cx2))/E;
          //nueva x
          std::complex<T> p = cx2 + h;

          if(std::abs(h) <= _precision){
            return p;
          }

          cx0 = cx1;
          cx1 = cx2;
          cx2 = p;

        }

        throw("Cantidad maxima de iteraciones alcanzada");

      }

  };

}

#endif
