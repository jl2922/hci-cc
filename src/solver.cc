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
#include "parallel.h"
#include "types.h"
#include "wavefunction.h"

namespace hci {

void Solver::load_wavefunction(const std::string& filename) {
  int n_in;
  double coef;
  const Det zero_det(n_orbs);
  int orb;
  
  // Read header line.
  std::ifstream wf_file(filename.c_str());
  wf_file >> n_in >> n_up >> n_dn;

  // Read each coef and det.
  for (int i = 0; i < n_in; i++) {
    wf_file >> coef;
    Det& det = wf.append_det(zero_det, coef);
    for (int j = 0; j < n_up; j++) {
      wf_file >> orb;
      det.up.set_orb(orb - 1, true);
    }
    for (int j = 0; j < n_dn; j++) {
      wf_file >> orb;
      det.dn.set_orb(orb - 1, true);
    }
  }

  Parallel::getInstance().checkpoint(
      "Loaded wavefunction with %d dets.\n", wf.size());
}

// Main solve procedure.
void Solver::solve() {
  setup();
  /// load_wavefunction("test/WAVE_small");
  load_wavefunction("test/WAVE_large");
  pt_det(0.00001);
}

// Deterministic 2nd-order purterbation.
void Solver::pt_det(const double eps_pt) {
  auto& parallel = Parallel::getInstance();

  // Save variational dets into hash set.
  std::unordered_set<Det, boost::hash<Det>> var_dets_set;
  const auto& var_dets = wf.get_dets();
  for (const auto& det: var_dets) {
    var_dets_set.insert(det);
  }
  var_dets_set.rehash(var_dets_set.size() * 10);

  // Estimate connections per det.
  const int n = wf.size();
  int hash_buckets = 0;
  const auto& var_coefs = wf.get_coefs();
  auto it_det = var_dets.begin();
  auto it_coef = var_coefs.begin();
  const int sample_interval = 100;
  for (int i = 0; i < n; i++) {
    const auto& det_i = *it_det++;
    const double coef_i = *it_coef++;
    if ((i % sample_interval) != 0) continue;
    const auto& connected_dets =
        find_connected_dets(det_i, eps_pt / fabs(coef_i));
    hash_buckets += connected_dets.size();
  }
  hash_buckets *= sample_interval / parallel.get_n_tasks() * 2;
  
  // Accumulate pt_sum for each det_a.
  std::unordered_map<Det, double, boost::hash<Det>> pt_sums;
  pt_sums.reserve(hash_buckets);
  it_det = var_dets.begin();
  it_coef = var_coefs.begin();
  int progress = 10;
  for (int i = 0; i < n; i++) {
    if (i % parallel.get_n_tasks() != parallel.get_task_id()) continue;
    const auto& det_i = *it_det++;
    const double coef_i = *it_coef++;
    const auto& connected_dets =
        find_connected_dets(det_i, eps_pt / fabs(coef_i));
    for (const auto& det_a: connected_dets) {
      if (var_dets_set.count(det_a) == 1) continue;
      const double H_ai = get_hamiltonian_elem(det_i, det_a, n_up, n_dn);
      if (fabs(H_ai) < Constants::EPSILON) continue;
      const double term = H_ai * coef_i;
      int task_owner_id =
          static_cast<int>(hash_value(det_a)) % parallel.get_n_tasks();
      if (task_owner_id == parallel.get_task_id()) {
        if (pt_sums.count(det_a) == 1) {
          pt_sums[det_a] += term;
        } else {
          pt_sums[det_a] = term;
        }
      } else {
        parallel.send<>(task_owner_id, DetDouble(det_a, term));
      }
    }

    if ((i + parallel.get_n_tasks()) * 100 >= n * progress) {
      // printf("Task %d: here prog %d\n", parallel.get_task_id(), progress);
      printf("Task %d Progress: %d%% (%d/%d) PT dets: %lu, hash load: %.2f\n",
          parallel.get_task_id(), progress, i, n, pt_sums.size(), pt_sums.load_factor());
      progress += 10;
    }
  }
  const auto& finishSign = DetDouble(Det(n_orbs), 0.0);
  parallel.broadcast<>(finishSign);
  int finishedTasks = 10;
  while(finishedTasks < parallel.get_n_tasks() - 1) {
    for (const auto& keyVal: parallel.collect<DetDouble>(finishSign)) {
      const auto& det_a = keyVal.first;
      if (det_a.is_zero()) finishedTasks++;
      const double term = keyVal.second;
      if (pt_sums.count(det_a) == 1) {
        pt_sums[det_a] += term;
      } else {
        pt_sums[det_a] = term;
      }
    }
  }
  parallel.checkpoint("All PT dets mapped to corresponding tasks.\n");
  printf("Task %d: Number of PT dets: %lu\n", parallel.get_task_id(), pt_sums.size());
  int n_pt_dets = parallel.sum(static_cast<int>(pt_sums.size()));
  parallel.checkpoint("Total Number of PT dets: %lu\n", n_pt_dets);

  // Accumulate contribution from each det_a to the pt_energy.
  double pt_energy = 0.0;
  for (auto it = pt_sums.begin(); it != pt_sums.end(); it++) {
    const auto& det_a = it->first;
    const double sum_a = it->second;
    const double H_aa = get_hamiltonian_elem(det_a, det_a, n_up, n_dn);
    pt_energy += pow(sum_a, 2) / (var_energy - H_aa);
  }
  printf("Task %d: PT correction: %.10f eV\n", parallel.get_task_id(), pt_energy);
  pt_energy = parallel.sum(pt_energy);
  parallel.checkpoint("Total PT correction: %.10f eV\n", pt_energy);
}
  
}
