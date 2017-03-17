#ifndef secant_f
#define secant_f

#include <cmath>
#include <limits>

template<class T>
class secant {


  private:
    //zero value
    T ZERO;
    //maximum number of iterations
    unsigned int MAXITER;
    //function
    T(*f_x)(T x);

    T slope(T x1, T x2, T y1, T y2){
      return (y2-y1)/(x2-x1);
    }

    T intersection(T x1, T y1, T m){
      return y1-(m*x1);
    }

    T distance(T x1, T x2, T y1, T y2){
      T x = x2-x1;
      T y = y2-y1;
      return (std::sqrt((x*x)+(y*y)));
    }


  public:

    secant(T(*pf_x)(T x), T zero = std::sqrt(std::numeric_limits<T>::epsilon()), int maxiter = 1000){
      ZERO = zero;
      MAXITER = maxiter;
      f_x = pf_x;
    }

    inline T operator() (T x1, T x2){
      //iteration counter
      unsigned int iterations = 0;
      //iteration variables
      //slope
      T m = 0;
      //intersection
      T b = 0;
      //X = 0
      T x0 = 0;
      //y = f_X(x0)
      T y0 = 0;

      //Y's calculation
      T y1 = f_x(x1), y2 = f_x(x2);

      //se valida que el menor este en x1
      if(y1 > y2){
        T tmp = x2;
        x2 = x1;
        x1 = tmp;
      }

      for(iterations; iterations < MAXITER; iterations++){
        //linear function calculos
        m = slope(x1, x2, y1, y2);
        b = intersection(x1,y1, m);
        //intersection with x
        x0 = -1*(b/m);
        //valor de y
        y0 = f_x(x0);

        //end determination
        if(fabs(y0) <= ZERO){
          return x0;
        }
        else{
          //se determina cual x se mantiene
          T d1 = distance(x1,y1,x0,y0);
          T d2 = distance(x2,y2,x0,y0);

          if(d1 < d2){
            x2 = x0;
            y2 = y0;
          }else{
            x1 = x0;
            y1 = y0;
          }
        }
      }

      return std::numeric_limits<T>::quiet_NaN();
    }

};


#endif
