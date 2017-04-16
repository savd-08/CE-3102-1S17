
#include "qr_decomp.hpp"
#include <limits>

using namespace anpi;




int main(int argc, char *argv[]){
	std::cout.precision(std::numeric_limits<double>::digits10 + 1);

	double d_arr[] = {1.0,2.0,2.0,1.0};
	Matrix<double> A(2, 2, d_arr);
    Matrix<double> Q;
	Matrix<double> R;
	std::vector<double> x = {0.0, 0.0};
	std::vector<double> b = {1.0, 1.0};
	qr(A, Q, R);

	print_Mat(A, 1);
	print_Mat(Q, 2);
	print_Mat(R, 3);
	double norm = testQR(A, Q, R);
	std::cout << "L2 norm" << norm << std::endl;


	solveQR(A, x, b);
	
	print_vec(x, 1);
    return 0;
}