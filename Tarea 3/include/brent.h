#ifndef secant_f
#define secant_f

#include <cmath>
#include <iostream>
#include <limits>

template<class T>
class brent {
  private:
    //maximum number of iterations
    unsigned int MAXITER;
    //zero value
    T ZERO;
    //function
    T(*f_x)(T x);

    void swap (T& t1, T& t2){
      T tmp = t1;
      t1 = t2;
      t2 = tmp;
    }

  public:
    brent(T(*pf_x)(T x), T pzero = std::sqrt(std::numeric_limits<T>::epsilon()), int maxiter = 1000){
        ZERO = pzero;
        MAXITER = maxiter;
        f_x = pf_x;
    }

    inline T operator()(T x1, T x2){
      //iterations counter
      unsigned int iterations = 0;
      //bandera de control
      bool ctrlflag = true;
      //calculo de las y
      T y1 = f_x(x1);
      T y2 = f_x(x2);
      //calculo de la tercera x
      T x3 = x1;
      //verifica que la menor este en x1
      if(fabs(y1) < fabs(y2)){
        swap(x1, x2);
        swap(y1, y2);
      }

      std::cout << "x1: " << x1 << " x2: " << x2;

      //calculo de la tercer punto
      x3 = x1;
      T y3 = y1;
      //raices actuales
      T x0 = 0;
      T y0 = 0;
      T x4 = 0;

      //iteraion
      for(iterations; iterations < MAXITER; iterations++){
        //si el error es menor que la tolerancia
        if(std::abs(x2-x1) < ZERO){

          return x0;
        }
        if(y1 != y3 && y2 != y3){
          //iterpolacion inversa cuadratica
          x0 =    (x1*y2*y3)/((y1-y2)*(y1-y3))
                + (x2*y1*y3)/((y2-y1)*(y2-y3))
                + (x3*y1*y2)/((y3-y1)*(y3-y2));
        }

        else{
          //secante
          x0 = (x2 - y2)*((x2 - x1)/(y2-y1));
        }
        //calculo de la y en la nueva raiz
        y0 = f_x(x0);
        //verificacion si hay que hacer biseccion
        if (    ( (x0 < (3 * x1 + x2) * 0.25) || (x0 > x2) ) ||
                ( ctrlflag && (std::abs(x0-x2) >= (std::abs(x2-x3) * 0.5)) ) ||
                ( !ctrlflag && (std::abs(x0-x2) >= (std::abs(x3-x4) * 0.5)) ) ||
                ( ctrlflag && (std::abs(x2-x3) < ZERO) ) ||
                ( !ctrlflag && (std::abs(x3-x4) < ZERO))  )
        {
            x0 = (x1+x2)*0.5;
            //cambio en la bandera
            ctrlflag = true;
        }
        else
        {
            ctrlflag = false;
        }

        //calculo de valores de final de iteraciones
        y0 = f_x(x0);

        x4 = x3;
        x3 = x2;
        //verifica cual si son de diferente signo
        if((y1*y0) < 0 ){
          x2 = x0;
        }
        else{
          x1 = x0;
        }
        //calculo de y's
        y1 = f_x(x1);
        y2 = f_x(x2);
        //verficia cual debe sera la primera
        if(fabs(y1) < fabs(y2)){
          swap(x1, x2);
          swap(y1, y2);
        }

      }
      return std::numeric_limits<T>::quiet_NaN();

    }

};

#endif
