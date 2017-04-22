#include "qr_decomp.hpp"
#include "lu_decomp.hpp"

int main(int argc, char *argv[]){
	std::cout.precision(std::numeric_limits<double>::digits10 + 1);

	double a1[] = {1,2,2,1};								//Sol: {1/3,1/3}
	double a2[] = {3.0,2.0,1.0,5.0,3.0,4.0,1.0,1.0,-1.0};	//Sol: {-4,6,1}
	double a3[] = {6,3,-1,1,3,1,1,1,1,2,1,0,2,3,4,0};		//Sol: {3,-4,1,2}
	double a4[] = {1,2,3,4,5,6,7,8,9}; 						//Singular
	double a5[] = {0.96,68,0.69,69,0.24,56,0.18,52,0.58,95,0.71,92,0.90,107,0.81,142};//Mal condicionada
	double a6[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

	anpi::Matrix<double> A(2, 2, a1); //Variar el vector de datos y las dimensiones
	std::vector<double> x(2); //Cambiar el tamaño de acuerdo con dimensiones de la matriz

	anpi::Matrix<double> A2(3, 3, a2);//Segunda Matriz de prueba
	std::vector<double> x2(3);

	//Vectores b para distintas matrices A
	std::vector<double> b1 = {1, 1};
	std::vector<double> b2 = {1.0,2.0,1.0};
	std::vector<double> b3 = {3,4,-4,-2};
	std::vector<double> b4 = {10,11,12};
	std::vector<double> b5 = {4,0.11,25,0.5};

	//------------------------------------------PRUEBAS QR------------------------------------------
	std::cout << "Descomposición QR" << std::endl;

	std::cout << "***************Prueba QR1***********************" << std::endl;

  anpi::Matrix<double> Q;
	anpi::Matrix<double> R;

	anpi::Matrix<double> Q2;
	anpi::Matrix<double> R2;

	//******Prueba de descomposición******
	anpi::qr(A, Q, R);

	anpi::print_Mat(A, 'A');
	std::cout << std::endl;
	anpi::print_Mat(Q, 'Q');
	std::cout << std::endl;
	anpi::print_Mat(R, 'R');
	std::cout << std::endl;

	//Comparación de Ax y Q*R
	double norm = anpi::testQR(A, Q, R);
	std::cout << "||A-Q*R||² = " << norm << std::endl;

	//*****Resolución de sistemas de ecuaciones*****
	std::cout << std::endl << "A*x = b" << std::endl << std::endl;
	anpi::print_Mat(A, 'A');
	std::cout << std::endl;
	anpi::print_vec(b1, 'b');//Cambiar ID de b
	std::cout << std::endl;

	anpi::solveQR(A, x, b1);//Cambiar ID de b
	anpi::print_vec(x, 'x');

	std::cout << "***************Prueba QR2***********************" << std::endl;

	//******Prueba de descomposición******
	anpi::qr(A2, Q2, R2);

	anpi::print_Mat(A2, 'A');
	std::cout << std::endl;
	anpi::print_Mat(Q2, 'Q');
	std::cout << std::endl;
	anpi::print_Mat(R2, 'R');
	std::cout << std::endl;

	//Comparación de Ax y Q*R
	double norm2 = anpi::testQR(A2, Q2, R2);
	std::cout << "||A-Q*R||² = " << norm2 << std::endl;

	//*****Resolución de sistemas de ecuaciones*****
	std::cout << std::endl << "A*x = b" << std::endl << std::endl;
	anpi::print_Mat(A2, 'A');
	std::cout << std::endl;
	anpi::print_vec(b2, 'b');//Cambiar ID de b
	std::cout << std::endl;

	anpi::solveQR(A2, x2, b2);//Cambiar ID de b
	anpi::print_vec(x2, 'x');





	//------------------------------------------PRUEBAS LU------------------------------------------

	std::cout << std::endl << "---------------------------------------------------------------------------------" << std::endl;
	std::cout << "Descomposición LU" << std::endl;

	std::cout << "***************Prueba LU1***********************" << std::endl;

	anpi::Matrix<double> A21(2, 2, a1);

	std::vector<double> x21(2);
	print_Mat(A21, 'A');
	/*******************************************************************************/

	//Inversion de matriz por descomposicion LU
	anpi::Matrix<double> Ai(2, 2, a6);
	anpi::invert(A21, Ai);
	std::cout << "Matriz A invertida: " << std::endl;
	anpi::print_Mat(Ai, 'i');

	/*******************************************************************************/

	//Solucion de Ax = b por descomposicion LU
	bool resolve = anpi::solveLU(A21, x21, b1, 2);
	if (resolve) {
			std::cout << "Con b = [1, 1]: " << std::endl;
			anpi::print_vec(x21, 'x');
	}
	else std::cout << "No se puede resolver el sistema lineal..." << std::endl;

	std::cout << "***************Prueba LU2***********************" << std::endl;

  anpi::Matrix<double> A22(3, 3, a2);

  std::vector<double> x22(3);
		print_Mat(A22, 'A');
  /*******************************************************************************/

  //Inversion de matriz por descomposicion LU
  anpi::Matrix<double> Ai2(3, 3, a6);
  anpi::invert(A22, Ai2);
  std::cout << "Matriz A invertida: " << std::endl;
  anpi::print_Mat(Ai2, 'i');

  /*******************************************************************************/

  //Solucion de Ax = b por descomposicion LU
  bool resolve2 = anpi::solveLU(A22, x22, b2, 3);
  if (resolve2) {
      std::cout << "Con b = [1, 2, 1]: " << std::endl;
      anpi::print_vec(x22, 'x');
  }
  else std::cout << "No se puede resolver el sistema lineal..." << std::endl;

	std::cout << std::endl << "---------------------------------------------------------------------------------" << std::endl;
	std::cout << "Matriz mal condicionada" << std::endl;

	anpi::Matrix<double> A3(4, 4, a5);//Segunda Matriz de prueba
	std::vector<double> x3(4);

	anpi::Matrix<double> Q3;
	anpi::Matrix<double> R3;

	//******Prueba de descomposición******
	anpi::qr(A3, Q3, R3);

	anpi::print_Mat(A3, 'A');
	std::cout << std::endl;
	anpi::print_Mat(Q3, 'Q');
	std::cout << std::endl;
	anpi::print_Mat(R3, 'R');
	std::cout << std::endl;

	//Comparación de Ax y Q*R
	double norm3 = anpi::testQR(A3, Q3, R3);
	std::cout << "||A-Q*R||² = " << norm << std::endl;

	//*****Resolución de sistemas de ecuaciones*****
	std::cout << std::endl << "A*x = b" << std::endl << std::endl;
	anpi::print_Mat(A3, 'A');
	std::cout << std::endl;
	anpi::print_vec(b5, 'b');//Cambiar ID de b
	std::cout << std::endl;

	anpi::solveQR(A3, x3, b5);//Cambiar ID de b
	anpi::print_vec(x3, 'x');

	std::cout << std::endl << "---------------------------------------------------------------------------------" << std::endl;
	std::cout << "Matriz con dimensiones incorrectas" << std::endl;

	anpi::Matrix<double> A4(8, 2, a5);
	std::vector<double> x4(8);

	anpi::Matrix<double> Q4;
	anpi::Matrix<double> R4; 

	try{
		anpi::qr(A4, Q4, R4);
	}
	catch(const std::exception& e){
		std::cerr << e.what() << std::endl;
	}

	std::cout << std::endl << "---------------------------------------------------------------------------------" << std::endl;
	std::cout << "Matriz singular" << std::endl;

	anpi::Matrix<double> A5(3, 3, a4);
	std::vector<double> x5(3);

	anpi::print_Mat(A5, 'A');

	try{
		anpi::solveQR(A5, x5, b4);
	}
	catch(const std::exception& e){
		std::cerr << e.what() << std::endl;
	}

  return 0;
}
