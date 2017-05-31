#ifndef _ANPI_LIEBMANN_HPP
#define _ANPI_LIEBMANN_HPP

#include <cmath>
#include <utility>
#include "Matrix.hpp"
#include <omp.h>
#include <iostream>
#include <fstream>

namespace anpi
{
  /**
   * Calcula la distribución de temperaturas.
   *
   * Parametros:
   * M matriz en la que se va a calcular la distribución. Se incializa en cero.
   * topT temperatura que va a estar constante en la parte superior de la matriz esta no se cuenta dentro de la matriz, para aislado debe ser un NaN.
   * rightT temperatura que va a estar constante en la parte derecha de la matriz esta no se cuenta dentro de la matriz, para aislado debe ser un NaN.
   * leftT temperatura que va a estar constante en la parte izquierda de la matriz esta no se cuenta dentro de la matriz, para aislado debe ser un NaN.
   * bottomT temperatura que va a estar constante en la parte inferior de la matriz esta no se cuenta dentro de la matriz, para aislado debe ser un NaN.
   * lambda valor para el metodo de sobre relajamiento.
   * tol porcentaje de error maximo que debe tener cada temperatura de la matriz
   */
  void liebmann(anpi::Matrix<double> &M, double topT, double rightT, double bottomT, double leftT, double lambda, double tol, bool enableOMP){


    M.fill(0.0);

    size_t lastRow = M.rows() - 1;
    size_t lastCol = M.cols() - 1;
    size_t tolCount = 0;
    size_t i;
    size_t j;

    double top  = 0.0, bottom = 0.0, left = 0.0, right = 0.0, newT = 0.0, oldT = 0.0, error = 0.0;

    while(tolCount < M.rows()*M.cols()){
      tolCount = 0;

      #pragma omp parallel for if(enableOMP) num_threads(2*omp_get_num_procs())\
      shared(M, lastRow, lastCol, tolCount, topT, rightT, bottomT, leftT, lambda, tol) \
      firstprivate(top, bottom, left, right, newT, oldT, error) \
      private(i,j) schedule(dynamic, static_cast<int>(M.rows()/omp_get_num_threads()))

      for (i = 0; i < M.rows(); i++) {
        for (j = 0; j < M.cols(); j++) {

          //condiciones para el de arriba y abajo
          if(i == 0){
            #pragma omp flush(M)
            bottom = M(i+1,j);
            if(topT != topT){
              top = bottom;
            }else{
              top = topT;
            }
          }else if(i == lastRow){
            #pragma omp flush(M)
            top = M(i-1,j);
            if(bottomT != bottomT){
              bottom = top;
            }else{
              bottom = bottomT;
            }
          }else{
            #pragma omp flush(M)
            top = M(i+1, j);
            bottom = M(i-1, j);
          }

          //condiciones para el derecho e izquierdo
          if(j == 0){
            #pragma omp flush(M)
            right = M(i,j+1);
            if(leftT != leftT){
              left = right;
            }else{
              left = leftT;
            }
          }else if(j == lastCol){
            #pragma omp flush(M)
            left = M(i,j-1);
            if(rightT != rightT){
              right = left;
            }else{
              right = rightT;
            }
          }else{
            #pragma omp flush(M)
            right = M(i, j+1);
            left = M(i, j-1);
          }

          //calculo de la nueva temperatura
          #pragma omp flush(M)
          oldT = M(i,j);
          newT = (top + right + bottom + left)/4.0;
          newT = (lambda*newT) + ((1.0-lambda)*oldT);
          #pragma omp atomic write
          M(i,j) = newT;

          //calculo de error
          error = std::abs((newT-oldT)/newT)*100;
          if(error <= tol){
            #pragma omp atomic
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
  void heatflux(anpi::Matrix<double> &M, anpi::Matrix< std::pair<double,double> > &M2, double topT, double rightT, double bottomT, double leftT, double k, bool enableOMP){
    M2.fill(std::pair<double,double>(0.0,0.0));

    size_t lastRow = M.rows() - 1;
    size_t lastCol = M.cols() - 1;
    size_t i;
    size_t j;

    double top  = 0.0, bottom = 0.0, left = 0.0, right = 0.0, qx = 0.0, qy = 0.0;

    #pragma omp parallel for if(enableOMP) num_threads(2*omp_get_num_procs())\
    shared(M, M2, lastRow, lastCol, topT, rightT, bottomT, leftT, k) \
    firstprivate(top, bottom, left, right, qx, qy) \
    private(i,j) schedule(dynamic, static_cast<int>(M.rows()/omp_get_num_threads()))

    for (i = 0; i < M.rows(); i++) {
      for (j = 0; j < M.cols(); j++) {
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
          top = M(i-1, j);
          bottom = M(i+1, j);
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
        //#pragma omp atomic write
        M2(i,j) = std::pair<double,double>(qx,qy);
      }//fin for j
    }//fin for i
  }

  void writeHeatMap(anpi::Matrix<double> &M, int vectorflag){
    std::ofstream heatmap;
    heatmap.open("ui/mapa_calor.txt");
    heatmap << vectorflag << "\n";

    for (size_t i = 0; i < M.rows(); i++) {
      for (size_t j = 0; j < M.cols(); j++) {
        heatmap << M(i,j) << "#";
      }
      heatmap << "\n";
    }
    heatmap.close();
  }

  void writeFlux(anpi::Matrix< std::pair<double,double> > M){
    std::ofstream qx;
    qx.open("ui/mapa_x.txt");

    std::ofstream qy;
    qy.open("ui/mapa_y.txt");

    for (size_t i = 0; i < M.rows(); i++) {
      for (size_t j = 0; j < M.cols(); j++) {
        qx << M(i,j).first << "#";
        qy << M(i,j).second << "#";
      }
      qx << "\n";
      qy << "\n";
    }
    qx.close();
    qy.close();
  }

}//fin de anpi

#endif
