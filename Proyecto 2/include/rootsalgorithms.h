#include <complex>
#include <limits>
#include <cmath>
#include "boost_poly.h"

namespace boost{ namespace math{ namespace tools{

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
      muller(polynomial<std::complex<T>> polinomio,T precision = std::sqrt(std::numeric_limits<T>::epsilon()), unsigned int iters = 1000){
        _precision = precision;
        _iters = iters;
        _polinomio = polinomio;
      }

      polynomial<std::complex<T>> operator()(T x1, T x2){

        //paso a numeros complejos
        std::complex<T> cx0 = x1;
        std::complex<T> cx2 = x2;

        std::complex<T> cy0 = poly_evaluator(cx0, _polinomio);//f_x(cx0);
        std::complex<T> cy2 = poly_evaluator(cx2, _polinomio);

        //secante para calcular un tercer punto
        std::complex<T> cx1 = (cx2 - cy2)*((cx2 - cx0)/(cy2 - cy0));

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

          if(std::abs(h) <= _precision){
            if(std::abs(std::real(p)) <= _precision){
              p = std::complex<T>(0.0, std::imag(p));
            }
            if(std::abs(std::imag(p)) <= _precision){
              p = std::complex<T>(std::real(p), 0.0);
            }

            p *= -1.0;
            polynomial<std::complex<T>> root{{p, std::complex<T>(1.0, 0.0)}};

            return root;
          }

          cx0 = cx1;
          cx1 = cx2;
          cx2 = p;

        }

        throw("Cantidad maxima de iteraciones alcanzada");

      }

  };

} //namespace tools
} //namespace math
} //namespace boost
