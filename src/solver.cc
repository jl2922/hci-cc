#include "solver.h"

#include <cmath>
#include <cstdio>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include "constants.h"
#include "det.h"
#include "wavefunction.h"

namespace hci {

// Main solve procedure.
void Solver::solve() {
  setup();
  wf.load("test/WAVE_small", n_orbs);
  pt_det(0.00001);
  printf("PT Energy: %.10f eV\n", pt_energy);
}

// Deterministic 2nd-order purterbation.
void Solver::pt_det(const double eps_pt) {
  std::unordered_set<hci::Det, boost::hash<hci::Det>> var_dets;
  std::list<hci::Det> pt_dets;
  std::unordered_map<hci::Det, double, boost::hash<hci::Det>> pt_sums;
  const int n = wf.n;

  // Save variational dets.
  for (int i = 0; i < n; i++) {
    const auto& det_i = wf.get_det(i);
    var_dets.insert(det_i);
  }

  // Calculate sums for each connected det_a.
  for (int i = 0; i < n; i++) {
    const auto& det_i = wf.get_det(i);
    const double coef_i = wf.get_coef(i);
    const auto& connected_dets =
        find_connected_dets(det_i, eps_pt / fabs(coef_i));
    // exit(0);
    for (int a = 0; a < connected_dets.n; a++) {
      const auto& det_a = connected_dets.get_det(a);
      if (var_dets.count(det_a) == 1) continue;
      const double H_ai = get_hamiltonian_elem(det_i, det_a);
      if (fabs(H_ai) < Constants::EPSILON) continue;
      const double term = H_ai * coef_i;
      if (pt_sums.count(det_a) == 1) {
        pt_sums[det_a] += term;
      } else {
        pt_dets.push_back(det_a);
        pt_sums[det_a] = term;
      }
    }
  }

  // Accumulate contribution from each det_a to the pt_energy.
  pt_energy = 0.0;
  for (auto it = pt_dets.begin(); it != pt_dets.end(); it++) {
    const auto& det_a = *it;
    const double sum_a = pt_sums[det_a];
    const double H_aa = get_hamiltonian_elem(det_a, det_a);
    pt_energy += pow(sum_a, 2) / (var_energy - H_aa);
  }
  printf("Number of PT dets: %d\n", pt_dets.size());
}
  
}
