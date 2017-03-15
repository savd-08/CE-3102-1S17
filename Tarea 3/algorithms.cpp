#include <cmath>
#include <iostream>

template <typename T>
T slope(T x1, T x2, T y1, T y2){
  return (y2-y1)/(x2-x1);
}

template <typename T>
T intersection(T x1, T y1, T m){
  return (m*x1)-y1;
}

template <typename T>
T distance(T x1, T x2, T y1, T y2){
  T x = x2-x1;
  T y = y2-y1;
  return (std::sqrt((x*x)+(y*y)));
}

//Secant Method
template <typename T>
T secant(T x1, T x2, T(*f_x)(T)){

  //zero value
  const T ZERO = 0.00001;
  //maximum number of iterations
  const unsigned int MAXITER = 50;
  //End flag
  bool end_flag = true;
  //iteration counter
  unsigned int i = 0;

  //iteration variables
  //slope
  T m = 0;
  //intersection
  T b = 0;
  //X = 0
  T x0 = 0;
  //y = f_X(x0)
  T y0 = 0;

  // se valida que se encierre una raiz
  if((x1*x2) > 0){
    throw "There is no root bracked";
  }
  //se valida que el menor este en x1
  if(x1 > x2){
    T tmp = x2;
    x2 = x1;
    x1 = tmp;
  }

  //Y's calculation
  T y1 = f_x(x1), y2 = f_x(x2);

  //loop
  while(i < MAXITER && end_flag){
    //linear function calculos
    m = slope(x1, x2, y1, y2);
    b = intersection(x1,y1, m);
    //intersection with x
    x0 = -1*(b/m);
    //valor de y
    y0 = f_x(x0);

    //end determination
    if(y0 <= ZERO || y0 >= (ZERO*-1)){
      end_flag = false;
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

    i++;
  }

  if(end_flag){
    throw "the iterations quantity has been exceeded";
  }

  return x0;

}

double function(double x){
  return (x*x + x);
}

int main(int argc, char const *argv[]) {

    double x1 = -0.75;
    double x2 = 2.39;

    double result = secant(x1, x2, &function);

    std::cout << "la raiz esta en la x: " << result << "\n";

    std::cout << "f(x): " << function(result) << "\n";


  return 0;
}
