#include <cmath>
#include <utility>
#include "Matrix.hpp"

/**
 * Calcula la distribución de temperaturas.
 *
 * Parametros:
 * M matriz en la que se va a calcular la distribución, esta va a ser llenada de ceros.
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

/**
*
*
*/
void heatflux(anpi::Matrix<double> &M, anpi::Matrix< std::pair<double,double> > &M2, double topT, double rightT, double bottomT, double leftT, double k){
  M2.fill(std::pair<double,double>(0.0,0.0));

  size_t lastRow = M.rows() - 1;
  size_t lastCol = M.cols() - 1;

  double top  = 0.0, bottom = 0.0, left = 0.0, right = 0.0, qx = 0.0, qy = 0.0, qn = 0.0, angle = 0.0;

  for (size_t i = 0; i < M.rows(); i++) {
    for (size_t j = 0; j < M.cols(); j++) {
      //condiciones para el de arriba y abajo
      if(i == 0){
        bottom = M(i+1,j);
        if(topT != topT){
          qy = 0;
        }else{
          qy = (-k)*((topT - bottom)/20);
        }
      }else if(i == lastRow){
        top = M(i-1,j);
        if(bottomT != bottomT){
          qy = 0;
        }else{
          qy = (-k)*((top - bottomT)/20);
        }
      }else{
        top = M(i+1, j);
        bottom = M(i-1, j);
        qy = (-k)*((top - bottom)/20);
      }
      //condiciones para el derecho e izquierdo
      if(j == 0){
        right = M(i,j+1);
        if(leftT != leftT){
          qx = 0;
        }else{
          qx = (-k)*((right-leftT)/20);
        }
      }else if(j == lastCol){
        left = M(i,j-1);
        if(rightT != rightT){
          qx = 0;
        }else{
          qx = (-k)*((rightT-left)/20);
        }
      }else{
        right = M(i, j+1);
        left = M(i, j-1);
        qx = (-k)*((rightT-leftT)/20);
      }
      //calculo de magnitud y angulo
      qn = std::sqrt(qx*qx + qy*qy);
      angle = std::atan2(qy,qx)*(180.0/3.1416);
      M2(i,j) = std::pair<double,double>(qn,angle);
    }//fin for j
  }//fin for i
}