
#include "qr_decomp.hpp"
#include "lu_decomp.hpp"
#include <limits>

using namespace anpi;	

int main(int argc, char *argv[]){
	std::cout.precision(std::numeric_limits<double>::digits10 + 1);

	/*double d_arr[] = {1.0,2.0,2.0,1.0};
	Matrix<double> A(2, 2, d_arr);
    Matrix<double> Q;
	Matrix<double> R;
	std::vector<double> x = {0.0, 0.0};
	std::vector<double> b = {1.0, 1.0};
	qr(A, Q, R);

	print_Mat(A, 1);
	print_Mat(Q, 2);
	print_Mat(R, 3);
	//double norm = testQR(A, Q, R);
	//std::cout << "L2 norm" << norm << std::endl;


	solveQR(A, x, b);	
	print_vec(x, 1);*/


	/*******************PRUEBA DESCOMPOSICION LU************************************/


/*	double test_mat[] = {3.0,-0.1,-0.2,
					0.1,7.0,-0.3,
					0.3,-0.2,10.0};    
	double test_mat[] = {1.0, -2.0, 3.0, 
					5.0, 8.0, -1.0,
					2.0, 1.0, 1.0};    
	double test_mat[]  = {1.0, 2.0, 3.0,
 					 4.0, 5.0, 6.0,
 					 7.0, 8.0, 10.0};   */
 	double test_mat[]  = {1.0, -2.0, 3.0,
 					 4.0, 0.0, 6.0,
 					 2.0, -1.0, 3.0};

    double init_mat[] = {0.0,0.0,0.0,
    	             0.0,0.0,0.0,
    	             0.0,0.0,0.0,
    	             0.0,0.0,0.0};
    int m = 3;
    int n = 3;

    Matrix<double> A(m, n, test_mat);
    Matrix<double> LU(m, n, init_mat);
   
    std::vector<double> x = {0.0, 0.0, 0.0};
    std::vector<double> b = {1.0, -2.0, -1.0};

    /*******************************************************************************/
    //Factorizacion LU por método de Crout y testLU
    //lu(A, LU);
 
    /*******************************************************************************/

    //Inversion de matriz por descomposicion LU
    Matrix<double> Ai(m, n, init_mat);

    invert(A, Ai);
    std::cout << "Matriz A invertida: " << std::endl;
    print_Mat(Ai, 'i');     

    /*******************************************************************************/

    //Solucion de Ax = b por LU 
    bool resolve = solveLU(A, x, b);
    if (resolve) {
        std::cout << "Con b = [3, 3, 4]: " << std::endl;   
        print_vec(x, 'x');  
    }
    else std::cout << "No se puede resolver el sistema lineal..." << std::endl;    

    return 0;
}
