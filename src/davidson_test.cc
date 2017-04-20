#include "davidson.h"

#include <cstdio>
#include <functional>

class TestHamiltonian {
 public:
  TestHamiltonian(int gamma) { this->gamma = gamma; }
  double get(int i, int j) {
    if (i == j) return -1.0 / (2 * i + 1);
    return -1.0 / gamma / (i + j + 1);
  }

 private:
  int gamma;
};

int main(int argc, char** argv) {
  int N = 10;
  TestHamiltonian hamiltonian(N);

  std::function<double(int, int)> getHamiltonian =
      std::bind(&TestHamiltonian::get, &hamiltonian, std::placeholders::_1, std::placeholders::_2);

  Davidson davidson(getHamiltonian, N);

  printf("Lowest Eigenvalue: %.10f\n", davidson.get_lowest_eigenvalue());
  return 0;
}