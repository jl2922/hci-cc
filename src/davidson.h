#ifndef DAVIDSON_H_
#define DAVIDSON_H_

#include <Eigen/Dense>
#include <functional>

class Davidson {
 public:
  Davidson(std::function<double(int, int)>& getHamiltonian, const int n) {
    this->getHamiltonian = &getHamiltonian;
    this->n = n;
  }

  double get_lowest_eigenvalue();

 private:
  std::function<double(int, int)>* getHamiltonian;
  int n;
};

double Davidson::get_lowest_eigenvalue() {
  Eigen::MatrixXd matrix(n, n);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      matrix(i, j) = (*getHamiltonian)(i, j);
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(matrix);
  const auto eigenvalues = eigenSolver.eigenvalues();
  
  double lowest_eigenvalue = eigenvalues[0];
  for (int i = 1; i < eigenvalues.size(); i++) {
    if (eigenvalues[i] < lowest_eigenvalue) {
      lowest_eigenvalue = eigenvalues[i];
    }
  }

  return lowest_eigenvalue;
}

#endif