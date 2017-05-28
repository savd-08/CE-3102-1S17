#include <iostream>
#include <cmath>
#include "Matrix.hpp"

/**
 * Calcula la propagación de temperaturas.
 *
 * Parametros:
 * M matriz en la que se va a calcular la propagación, esta va a ser llenada de ceros.
 * topT temperatura que va a estar constante en la parte superior de la matriz esta no se cuenta dentro de la matriz, para aislado debe ser un NaN.
 * rightT temperatura que va a estar constante en la parte derecha de la matriz esta no se cuenta dentro de la matriz, para aislado debe ser un NaN.
 * leftT temperatura que va a estar constante en la parte izquierda de la matriz esta no se cuenta dentro de la matriz, para aislado debe ser un NaN.
 * bottomT temperatura que va a estar constante en la parte inferior de la matriz esta no se cuenta dentro de la matriz, para aislado debe ser un NaN.
 * lambda valor para el metodo de sobre relajamiento.
 * tol porcentaje de error maximo que debe tener cada temperatura de la matriz
 */
void liebmann(anpi::Matrix<double> &M, double topT, double rightT, double bottomT, double leftT, double lambda, double tol){

  M.fill(0.0);

  size_t lastRow = M.rows() - 1;
  size_t lastCol = M.cols() - 1;
  size_t tolCount = 0;

  double top  = 0.0, bottom = 0.0, left = 0.0, right = 0.0, newT = 0.0, oldT = 0.0, error = 0.0;

  while(tolCount < M.rows()*M.cols()){
    tolCount = 0;

    for (size_t i = 0; i < M.rows(); i++) {
      for (size_t j = 0; j < M.cols(); j++) {
        oldT = M(i,j);
        //condiciones para el de arriba y abajo
        if(i == 0){
          bottom = M(i+1,j);
          if(topT != topT){
            top = bottom;
          }else{
            top = topT;
          }
        }else if(i == lastRow){
          top = M(i-1,j);
          if(bottomT != bottomT){
            bottom = top;
          }else{
            bottom = bottomT;
          }
        }else{
          top = M(i+1, j);
          bottom = M(i-1, j);
        }
        //condiciones para el derecho e izquierdo
        if(j == 0){
          right = M(i,j+1);
          if(leftT != leftT){
            left = right;
          }else{
            left = leftT;
          }
        }else if(j == lastCol){
          left = M(i,j-1);
          if(rightT != rightT){
            right = left;
          }else{
            right = rightT;
          }
        }else{
          right = M(i, j+1);
          left = M(i, j-1);
        }
        //calculo de la nueva temperatura
        newT = (top+right+bottom+left)/4.0;
        newT = (lambda*newT) + ((1.0-lambda)*oldT);
        M(i,j) = newT;
        //calculo de error
        error = std::abs((newT-oldT)/newT)*100;
        if(error <= tol){
          tolCount++;
        }
      }//fin for j
    }//fin for i
  }//fin while
}

int main(int argc, char const *argv[]) {

  anpi::Matrix<double> M(3,3,0.0);

  liebmann(M, 100.0, 50.0, (0.0/0.0), 75.0, 1.5, 1.0);

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      std::cout << M(i,j) << " ";
    }
    std::cout << "\n";
  }

  return 0;
}
