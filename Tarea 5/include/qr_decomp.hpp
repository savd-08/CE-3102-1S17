#ifndef _ANPI_QR_DECOMP_H
#define _ANPI_QR_DECOMP_H

#include "linear_utils.hpp"
#include <algorithm> 
#include <math.h> 

namespace anpi{

	/*Descomposición QR a través de transformaciones de Householder,
	 *para una matriz cuadrada A, con punteros a Q y R como parámetros.
	 */
	template<typename T>
	bool qr(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R){
		int i,j,k;
		int n = A.rows();
		Q = Matrix<double>(n, n, T(0));
		R = A;
		bool sing = false;
		std::vector<T>  c(n), d(n);
		T scale, sigma, sum, tau;

		for (k = 0; k < n-1; k++) {
			scale = T(0);
			for (i = k; i < n; i++) {
				scale = std::max(scale, fabs(R[i][k]));
			}
			if (fabs(scale) <= std::sqrt(std::numeric_limits<T>::epsilon())) {
				sing = true;
				c[k] = d[k] = T(0);
			} else {


				//Cálculo de la norma para la columna Ri
				for (sum = T(0), i = k; i < n; i++) {
					R[i][k] /= scale;
					sum += R[i][k]*R[i][k];
				}
				sigma = std::sqrt(sum);


				sigma = (R[k][k] < 0) ? -sigma : sigma;


				R[k][k] += sigma;
				c[k] = sigma*R[k][k];
				d[k] = -scale*sigma;

				for (j = k + 1; j < n; j++) {

					for (sum = T(0), i = k; i < n; i++) {
						sum += R[i][k]*R[i][j];
					}

					tau = sum/c[k];
					for (i = k; i < n; i++) {
						R[i][j] -= tau*R[i][k];
					}
				}
			}
		}

		d[n - 1] = R[n - 1][n - 1];
		if (fabs(d[n-1]) <= std::sqrt(std::numeric_limits<T>::epsilon())) {
			sing = true;
		}

		//Se genera la matriz Qt (la misma que Q por ortogonalidad)
		for (i = 0; i < n; i++) {
			Q[i][i] = T(1);
		}

		for (k = 0; k < n - 1; k++) {
			if (c[k] != T(0)) {
				for (j = 0; j < n; j++) {
					sum = T(0);
					for (i = k; i < n; i++)
					sum += R[i][k]*Q[i][j];
					sum /= c[k];
					for (i = k; i < n; i++)
					Q[i][j] -= sum*R[i][k];
				}
			}
		}

		//Se genera la matriz R
		for (i = 0; i < n; i++) {

			R[i][i] = d[i];
			for (j = 0; j < i; j++) {
				R[i][j] = T(0);
			}
		}

		return sing;
	}

	template<typename T>
	T testQR(const Matrix<T>& A, Matrix<T>& Q, Matrix<T>& R){
		Matrix<T> A_x = Q*R;
		//Comparación de matrices a través de norma L2
		return matrix_compare(A, A_x); 
	}

	/*Soluciona un sistema de ecuaciones a través de descomposición QR
	 *utilizando transformaciones de Householder.
	 */
	template<typename T>
	void solveQR(const Matrix<T>& A, std::vector<T>& x, std::vector<T>& b) {
		Matrix<T> Q;
		Matrix<T> R;
		bool sing = qr(A, Q, R); //Descomposición QR
		b = Q*b;
		bckwd_subst(R, x, b, sing); //Sustitución hacia atrás para encontrar x
	}
} // namespace ANPI

#endif
