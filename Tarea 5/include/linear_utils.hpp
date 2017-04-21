#ifndef _ANPI_LINEAR_UTILS_H
#define _ANPI_LINEAR_UTILS_H

#include <limits>
#include "Matrix.hpp"
#include "Exception.hpp"
#include <vector>
#include <iostream>

namespace anpi{

	//Imprime matrices
	template<typename T>
	void print_Mat(const Matrix<T>& M, char id){
		int n = M.rows();
		int i, j;
		std::cout << "Mat " << id << std::endl;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				std::cout << M[i][j] << '\t';
			}
			std::cout << std::endl;
		}
	}

	//Imprime vectores
	template<typename T>
	void print_vec(const std::vector<T> v, char id){
		int n = v.size();
		int i;
		std::cout << "Vec " << id << std::endl;
		for (i = 0; i < n; i++) {
			std::cout << v[i] << '\t';
		}
		std::cout << std::endl;
	}

	//Norma L2 (||A-B||^2)
	template<typename T>
	T matrix_compare(const Matrix<T>& A, const Matrix<T>& B) {
	  if (A.rows() != B.rows() or  A.cols() != B.cols()){
	  	throw Exception("Las matrices no son de dimensiones iguales");
	  }
	  T res = 0;
	  for(int i = 0; i < A.rows(); i++) {
	    for (int j = 0; j < A.cols(); j++) {
	      res += (A[i][j]-B[i][j]) * (A[i][j]-B[i][j]);
	    }
	  }

	  res /= A.rows()*A.cols();
	  return res;
	}

	//Sustitución hacia atrás para encontrar x (solución de Rx=b)
	template<typename T>
	void bckwd_subst(const Matrix<T>& R, std::vector<T> &x, const std::vector<T> &b, bool sing) {
		int n = R.rows()-1;
		if (sing) throw Exception("La matriz A es singular");
		//Método de sustitución hacia atrás (ver notas del curso)
		x[n] = b[n]/R[n][n];
		for (int i = n - 1; i >= 0; i--)
		{
			x[i] = (b[i] - R[i][i + 1]*x[i + 1])/R[i][i];
		}
	}
} // namespace ANPI

#endif
