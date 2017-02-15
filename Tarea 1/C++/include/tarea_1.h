#include <iostream>
#include <math.h> 

using namespace std;

//Estructura utilizada para contener las soluciones de una ecuación cuadrática
template <typename fp> 
struct quad_ec_sol { 
  fp x1, x2; 
}; 

template <typename fp>
fp ej_2 (int, fp); 

template <typename fp>
quad_ec_sol<fp> ej_3 (fp, fp, fp);

