#include "solver.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <list>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include "constants.h"
#include "det.h"
#include "wavefunction.h"

namespace hci {

void Solver::load_wavefunction(const std::string& filename) {
  int n_in;
  double coef;
  const Det det_proto(n_orbs);
  int orb;
  
  // Read header line.
  std::ifstream wf_file(filename.c_str());
  wf_file >> n_in >> n_up >> n_dn;

  // Read each coef and det.
  for (int i = 0; i < n_in; i++) {
    wf_file >> coef;
    Det& det = wf.append_det(det_proto, coef);
    for (int j = 0; j < n_up; j++) {
      wf_file >> orb;
      det.up.set_orb(orb - 1, true);
    }
    for (int j = 0; j < n_dn; j++) {
      wf_file >> orb;
      det.dn.set_orb(orb - 1, true);
    }
  }

  printf("Loaded wavefunction with %d dets.\n", wf.size());
}

// Main solve procedure.
void Solver::solve() {
  setup();
  load_wavefunction("test/WAVE_small");
  // load_wavefunction("test/WAVE_large");
  pt_det(0.00001);
  printf("PT Energy: %.10f eV\n", pt_energy);
}

// Deterministic 2nd-order purterbation.
void Solver::pt_det(const double eps_pt) {
  std::unordered_set<hci::Det, boost::hash<hci::Det>> var_dets_set;
  std::list<hci::Det> pt_dets;
  std::unordered_map<hci::Det, double, boost::hash<hci::Det>> pt_sums;
  const int n = wf.size();

  // Save variational dets.
  const auto& var_dets = wf.get_dets();
  for (const auto& det: var_dets) {
    var_dets_set.insert(det);
  }

  // Calculate sums for each connected det_a.
  const auto& var_coefs = wf.get_coefs();
  auto it_det = var_dets.begin();
  auto it_coef = var_coefs.begin();
  int progress = 10;
  for (int i = 0; i < n; i++) {
    const auto& det_i = *it_det++;
    const double coef_i = *it_coef++;
    const auto& connected_dets =
        find_connected_dets(det_i, eps_pt / fabs(coef_i));
    for (const auto& det_a: connected_dets.get_dets()) {
      if (var_dets_set.count(det_a) == 1) continue;
      const double H_ai = get_hamiltonian_elem(det_i, det_a, n_up, n_dn);
      if (fabs(H_ai) < Constants::EPSILON) continue;
      const double term = H_ai * coef_i;
      if (pt_sums.count(det_a) == 1) {
        pt_sums[det_a] += term;
      } else {
        pt_dets.push_back(det_a);
        pt_sums[det_a] = term;
      }
    }
    if ((i + 1) * 100 >= n * progress) {
      printf("Progress: %d%% (%d/%d)\n", progress, i, n);
      progress += 10;
    }
  }

  // Accumulate contribution from each det_a to the pt_energy.
  pt_energy = 0.0;
  for (auto it = pt_dets.begin(); it != pt_dets.end(); it++) {
    const auto& det_a = *it;
    const double sum_a = pt_sums[det_a];
    const double H_aa = get_hamiltonian_elem(det_a, det_a, n_up, n_dn);
    pt_energy += pow(sum_a, 2) / (var_energy - H_aa);
  }
  printf("Number of PT dets: %d\n", static_cast<int>(pt_dets.size()));
}
  
}
