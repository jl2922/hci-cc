#include "davidson.h"

#include <cstdio>
#include <functional>
#include <list>

// Test with a Hilbert matrix.
class TestHamiltonian {
 public:
  TestHamiltonian(int n, int gamma) {
    this->n = n;
    this->gamma = gamma;
  }

  double get_hamiltonian(int i, int j) {
    if (i == j) return -1.0 / (2 * i + 1);
    return -1.0 / gamma / (i + j + 1);
  }

  std::vector<double> apply_hamiltonian(const std::vector<double>& v) {
    std::vector<double> Hv(n, 0.0);
    for (int i = 0; i < n; i++) {
      Hv[i] += get_hamiltonian(i, i) * v[i];
      for (int j = i + 1; j < n; j++) {
        double h_ij = get_hamiltonian(i, j);
        Hv[i] += h_ij * v[j];
        Hv[j] += h_ij * v[i];
      }
    }
    return Hv;
  }

 private:
  int n;
  int gamma;
};

int main(int argc, char** argv) {
  int N = 10000;
  TestHamiltonian hamiltonian(N, 10);

  std::function<double(int, int)> get_hamiltonian = std::bind(
      &TestHamiltonian::get_hamiltonian,
      &hamiltonian,
      std::placeholders::_1,
      std::placeholders::_2);

  std::function<std::vector<double>(std::vector<double>)> apply_hamiltonian =
      std::bind(&TestHamiltonian::apply_hamiltonian, &hamiltonian, std::placeholders::_1);

  Davidson davidson(get_hamiltonian, apply_hamiltonian, N);

  printf("Lowest Eigenvalue: %.10f\n", davidson.get_lowest_eigenvalue());

  const auto& eigenvector = davidson.get_lowest_eigenvector();

  printf("Top elements in eigenvector:\n");
  for (int i = 0; i < 5; i++) printf("%d: %.10f\n", i, eigenvector[i]);

  return 0;
}