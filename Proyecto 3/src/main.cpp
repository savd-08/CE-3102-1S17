#include <iostream>
#include "liebmann.hpp"



int main(int argc, char const *argv[]) {

  anpi::Matrix<double> M(3,3,0.0);

  liebmann(M, 100.0, 50.0, (0.0/0.0), 75.0, 1.5, 1.0);

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      std::cout << M(i,j) << " ";
    }
    std::cout << "\n";
  }

  std::cout << "\n\n\n";

  anpi::Matrix< std::pair<double,double> > M2(3,3,std::pair<double,double>(0.0,0.0));

  heatflux(M, M2, 100.0, 50.0, (0.0/0.0), 75.0, 0.49);

  for (size_t i = 0; i < 3; i++) {
    for (size_t j = 0; j < 3; j++) {
      std::cout << "(" << M2(i,j).first << "," << M2(i,j).second << ")" << " ";
    }
    std::cout << "\n";
  }

  return 0;
}
