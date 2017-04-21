#ifndef DAVIDSON_H_
#define DAVIDSON_H_

#include <Eigen/Dense>
#include <functional>
#include <vector>
#include <iostream>

// Translated from Adam's fortran code.
class Davidson {
 public:
  Davidson(
      std::function<double(int, int)>& get_hamiltonian_elem,
      std::function<std::vector<double>(std::vector<double>)>& apply_hamiltonian,
      const int n) {
    this->get_hamiltonian_elem = &get_hamiltonian_elem;
    this->apply_hamiltonian = &apply_hamiltonian;
    this->n = n;
    diagonalized = false;
  }

  void set_initial_vector(const std::vector<double>& vector) { this->initial_vector = vector; }

  void diagonalize();

  double get_lowest_eigenvalue() {
    if (!diagonalized) diagonalize();
    return lowest_eigenvalue;
  }

  const std::vector<double>& get_lowest_eigenvector() {
    if (!diagonalized) diagonalize();
    return lowest_eigenvector;
  }

 private:
  // Use functional programming to allow either direct or indirect evaluation.
  std::function<double(int, int)>* get_hamiltonian_elem;
  std::function<std::vector<double>(std::vector<double>)>* apply_hamiltonian;

  // Solutions.
  double lowest_eigenvalue;
  std::vector<double> lowest_eigenvector;

  std::vector<double> initial_vector;

  // Length in each direction.
  int n;

  bool diagonalized;
};

void Davidson::diagonalize() {
  if (n == 1) {
    lowest_eigenvalue = (*get_hamiltonian_elem)(0, 0);
    lowest_eigenvector = std::vector<double>(1, 1.0);
  }

  const int iterations = std::min(n, 50);
  double lowest_eigenvalue = 0.0;
  double lowest_eigenvalue_prev = 0.0;
  double residual_norm = 0.0;

  Eigen::MatrixXd v = Eigen::MatrixXd::Zero(n, iterations);

  if (static_cast<int>(initial_vector.size()) != n) {
    v(0, 0) = 1.0;  // Start from HF.
  } else {
    for (int i = 0; i < n; i++) v(i, 0) = initial_vector[i];
    v.col(0).normalize();
  }

  Eigen::MatrixXd Hv = Eigen::MatrixXd::Zero(n, iterations);
  Eigen::VectorXd w = Eigen::VectorXd::Zero(n);  // Lowest eigenvector so far.
  Eigen::VectorXd Hw = Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd h_krylov = Eigen::MatrixXd::Zero(iterations, iterations);
  Eigen::MatrixXd h_overwrite = Eigen::MatrixXd::Zero(iterations, iterations);
  Eigen::VectorXd eigenvalues = Eigen::VectorXd::Zero(iterations);
  int len_work = 3 * iterations - 1;
  Eigen::VectorXd work(len_work);
  bool converged = false;
  std::vector<double> tmp_v(n);

  // Get diagonal elements.
  Eigen::VectorXd diag_elems(n);
  for (int i = 0; i < n; i++) diag_elems[i] = (*get_hamiltonian_elem)(i, i);

  // First iteration.
  for (int i = 0; i < n; i++) tmp_v[i] = v(i, 0);
  const auto& tmp_Hv = (*apply_hamiltonian)(tmp_v);
  for (int i = 0; i < n; i++) Hv(i, 0) = tmp_Hv[i];
  lowest_eigenvalue = v.col(0).dot(Hv.col(0));
  h_krylov(0, 0) = lowest_eigenvalue;
  w = v.col(0);
  Hw = Hv.col(0);
  printf("Initial iteration. Eigenvalue: %.10f\n", lowest_eigenvalue);

  residual_norm = 1.0;  // So at least one iteration is done.
  int n_iter = std::min(n, iterations + 1);
  int n_diagonalize = 1;

  for (int it = 1; it < n_iter; it++) {
    // Compute residual.
    for (int j = 0; j < n; j++) {
      v(j, it) = (Hw(j, 0) - lowest_eigenvalue * w(j, 0)) / (lowest_eigenvalue - diag_elems(j));
      if (fabs(lowest_eigenvalue - diag_elems[j]) < 1.0e-8) v(j, it) = -1.0;
    }

    // If residual is small, converge.
    residual_norm = v.col(it).norm();
    if (residual_norm < 1.0e-12) converged = true;

    // Orthogonalize and normalize.
    for (int i = 0; i < it; i++) {
      double norm = v.col(it).dot(v.col(i));
      v.col(it) -= norm * v.col(i);
    }
    v.col(it).normalize();

    // Apply H once.
    for (int i = 0; i < n; i++) tmp_v[i] = v(i, it);
    const auto& tmp_Hv2 = (*apply_hamiltonian)(tmp_v);
    for (int i = 0; i < n; i++) Hv(i, it) = tmp_Hv2[i];

    // Construct Krylow matrix and diagonalize.
    for (int i = 0; i <= it; i++) {
      h_krylov(i, it) = v.col(i).dot(Hv.col(it));
      h_krylov(it, i) = h_krylov(i, it);
    }
    std::cout << h_krylov.leftCols(it).topRows(it) << std::endl;
    len_work = 3 * it + 2;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(h_krylov);
    const auto& eigenvalues = eigenSolver.eigenvalues();
    const auto& eigenvectors = eigenSolver.eigenvectors();
    lowest_eigenvalue = eigenvalues[0];
    int lowest_id = 0;
    for (int i = 1; i <= it; i++) {
      if (eigenvalues[i] < lowest_eigenvalue) {
        lowest_eigenvalue = eigenvalues[i];
        lowest_id = i;
      }
    }
    w = v.leftCols(it) * eigenvectors.col(lowest_id).topRows(it);
    Hw = Hv.leftCols(it) * eigenvectors.col(lowest_id).topRows(it);

    if (it > 1 && fabs(lowest_eigenvalue - lowest_eigenvalue_prev) < 1.0e-6) {
      converged = true;
      break;
    } else {
      lowest_eigenvalue_prev = lowest_eigenvalue;
      n_diagonalize++;
      printf("Iteration #%d. Eigenvalue: %.10f\n", n_diagonalize, lowest_eigenvalue);
    }
    if (converged) break;
  }

  this->lowest_eigenvalue = lowest_eigenvalue;
  lowest_eigenvector.resize(n);
  for (int i = 0; i < n; i++) lowest_eigenvector[i] = w(i);
  diagonalized = true;
}

#endif
