#ifndef _ANPI_LU_DECOMP_H
#define _ANPI_LU_DECOMP_H

#include "linear_utils.hpp"
#include <algorithm> 
#include <math.h> 
#include <cmath>

namespace anpi{

	/*Descomposición LU a través del método de Crout,
	 *para una matriz cuadrada A, y dos matrices L y U.
	 */
	template<typename T>
	void lu(const Matrix<T>& A, Matrix<T>& LU){

        T res = T(0);
        int n = A.rows();

        //Descomposicion LU por Crout
        for (int i = 0; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                res = T(0);
                for (int k = 0; k < i; k++)
                    res += LU[i][k] * LU[k][j];
                LU[i][j] = A[i][j] - res;
            }

            for (int j = i + 1; j < n; j++)
            {
                res = T(0);
                for (int k = 0; k < i; k++)
                    res += LU[j][k] * LU[k][i];
                LU[j][i] = (T(1) / LU[i][i]) * (A[j][i] - res);
            }
        }

        //Se desglosa LU en L y U 
        Matrix<T> L = LU;
        Matrix<T> U = LU;

        for ( int i = 0; i < n; i++) {
            for ( int j = 0; j < n; j++) {
               if (i == j) U[i][j] = T(1);
               if (i < j) L[i][j] = T(0);
               if (i > j) U[i][j] = T(0);
            }
        }

        print_Mat(A, 'A');
        print_Mat(L, 'L');
        print_Mat(U, 'U');
        
        testLU(A,L,U);
        
	}

    /**************************************************************************************/

    /* Test para calcular la diferencia de la norma entre la matriz reconstruida LU
    * y la matriz propuesta A para verificar cuántitativamente el valor en que difieren
    */
    template<typename T>
    T calc_norm(const Matrix<T>& mat){
        int n = mat.rows();
        T res = T(0);

        for ( int i = 0; i < n; i++) {
            for ( int j = 0; j < n; j++) {
               res = res + (std::pow(mat[i][j],2));
            }
        }
        return std::sqrt(res);
    }

	
	template<typename T>
	void testLU(const Matrix<T>& A, Matrix<T>& L, Matrix<T>& U){

        Matrix<T> LU = L * U; //Se construye la matriz LU
        T norm_LU = calc_norm(LU); //Se calcula la norma de LU
        T norm_A = calc_norm(A); //Se calcula la norma de A

        std::cout <<  "La diferencia entre la norma de la matriz reconstruida LU y A es: " <<
                       std::abs(norm_A - norm_LU) << std::endl;
        
	}


    /**************************************************************************************/

    /*Calcula la solucion de Ax = b 
     * recibiendo una matriz A y dos vectores x y b como parámetros
     */
	template<typename T>
	bool solveLU(const Matrix<T>& A, std::vector<T>& x, const std::vector<T>& b) {

            int n = A.rows();
            bool resolve = true;
            T arr[n];

            for (int i = 0; i < n; i++) {
                arr[i] = T(0);
            }

            Matrix<T> LU(n, n, arr);
    
            T sum = T(0);
            
            //Se descompone la matriz mediante LU
            lu(A, LU);    
            
            // soluciona Ly = b por medio de LU
            std::vector<T> y(n);
            for (int i = 0; i < n; i++)
            {
                sum = T(0);
                for (int k = 0; k < i; k++)
                    sum += LU[i][k] * y[k];
                y[i] = b[i] - sum;
                
            }

            // soluciona Ux = y por medio de LU
            for (int i = n - 1; i >= 0; i--)
            {
                sum = T(0);
                for (int k = i + 1; k < n; k++)
                    sum += LU[i][k] * x[k];
                x[i] = (T(1) / LU[i][i]) * (y[i] - sum);
                if (isnan(y[i])) {
                    resolve = false;
                    break;
                }
            }
            
            //indica si el sistema se puede resolver 
            return resolve;
    }

     

    /**************************************************************************************/

    
    
	template<typename T>
	void invert(const Matrix<T>& A, Matrix<T>& Ai) {

        int n = A.rows();
        std::vector<T> b(n), x(n);
        int i, j;

        for (int k = 0; k < n; k++) {
            x[k] = T(0);
        }
         
        for(i = 0; i < n; i++) {

            // Vector de términos independientes, columna i de la identidad
            for(j = 0; j < n; j++) b[j] = (j == i) ? 1 : 0;
               
            // Resuelve el problema y obtiene la columna i de la inversa
            bool sol = solveLU(A, x, b);
             
            // Guarda la columna x
            for(j = 0; j < n; j++) Ai[j][i] = x[j];
        }
        
    }   
     

} // namespace ANPI

#endif