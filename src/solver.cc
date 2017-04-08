#include "solver.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <mutex>
#include <iostream>
#include <list>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <boost/mpi.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/serialization/utility.hpp>
#include "big_unordered_map.h"
#include "constants.h"
#include "det.h"
#include "parallel.h"
#include "types.h"
#include "wavefunction.h"

namespace hci {

Solver::Solver() {
  mpi.id = mpi.world.rank();
  mpi.n = mpi.world.size();
}

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

  mpi.world.barrier();
  if (mpi.id == 0) printf("Loaded wavefunction w/ %d dets.\n", wf.size());
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
  // Save variational dets into hash set.
  std::unordered_set<Det, boost::hash<Det>> var_dets_set;
  const auto& var_dets = wf.get_dets();
  for (const auto& det: var_dets) {
    var_dets_set.insert(det);
  }
  var_dets_set.rehash(var_dets_set.size() * 10);

  // Estimate connections per det.
  const int n = wf.size();
  int hash_buckets_local = 0, hash_buckets = 0;
  const auto& var_coefs = wf.get_coefs();
  auto it_det = var_dets.begin();
  auto it_coef = var_coefs.begin();
  const int sample_interval = 100;
  for (int i = 0; i < n; i++) {
    const auto& det_i = *it_det++;
    const double coef_i = *it_coef++;
    if ((i % (sample_interval * mpi.n)) != sample_interval * mpi.id) continue;
    const auto& connected_dets =
        find_connected_dets(det_i, eps_pt / fabs(coef_i));
    hash_buckets_local += connected_dets.size();
  }
  hash_buckets_local *= 2 * sample_interval;
  all_reduce(mpi.world, hash_buckets_local, hash_buckets, std::plus<int>());
  
  // Accumulate pt_sum for each det_a.
  std::unordered_map<Det, double, boost::hash<Det>> pt_sums;
  std::pair<std::vector<BitsBlock>, double> skeleton;
  skeleton.first = Det(n_orbs).as_vector();
  BigUnorderedMap<std::vector<BitsBlock>, double, 
      boost::hash<std::vector<BitsBlock>>> pt_sums_new(mpi.world, skeleton, 200);
  pt_sums_new.reserve(hash_buckets);
  if (mpi.id == 0) printf("Setup hash map with # buckets: %d\n", hash_buckets);
  // pt_sums.reserve(hash_buckets);
  pt_sums_new.reserve(hash_buckets);
  it_det = var_dets.begin();
  it_coef = var_coefs.begin();
  int progress = 10;
  for (int i = 0; i < n; i++) {
    const auto& det_i = *it_det++;
    const double coef_i = *it_coef++;
    if (i % mpi.n != mpi.id) continue;
    const auto& connected_dets =
        find_connected_dets(det_i, eps_pt / fabs(coef_i));
    for (const auto& det_a: connected_dets) {
      if (var_dets_set.count(det_a) == 1) continue;
      const double H_ai = get_hamiltonian_elem(det_i, det_a, n_up, n_dn);
      if (fabs(H_ai) < Constants::EPSILON) continue;
      const double term = H_ai * coef_i;
      // pt_sums[det_a] += term;
      pt_sums_new.async_inc(det_a.as_vector(), term);
    }
    if ((i + 1) * 100 >= n * progress && mpi.id == 0) {
      const auto& local_map = pt_sums_new.get_local_map();
      printf("MASTER: Progress: %d%% (%d/%d), PT dets: %lu, hash load: %.2f\n",
          progress, i, n, local_map.size(), local_map.load_factor());
      progress += 10;
    }
  }
  pt_sums_new.complete_async_incs();

  std::size_t n_pt_dets = pt_sums_new.size();
  if (mpi.id == 0) printf("Number of PT dets: %lu\n", n_pt_dets);

  // printf("Number of PT dets: %lu\n", pt_sums.size());

  // Accumulate contribution from each det_a to the pt_energy.
  pt_energy = 0.0;
  Det det_a;
  double pt_energy_local = 0.0;
  for (const auto& kv: pt_sums_new.get_local_map()) {
    det_a.from_vector(kv.first, n_orbs);
    const double sum_a = kv.second;
    // printf("sum_a: %.10f\n", sum_a);
    const double H_aa = get_hamiltonian_elem(det_a, det_a, n_up, n_dn);
    // printf("H_aa: %.10f\n", H_aa);
    pt_energy_local += pow(sum_a, 2) / (var_energy - H_aa);
  }
  reduce(mpi.world, pt_energy_local, pt_energy, std::plus<double>(), 0);
  // for (auto it = pt_sums.begin(); it != pt_sums.end(); it++) {
  //   const auto& det_a = it->first;
  //   const double sum_a = it->second;
  //   const double H_aa = get_hamiltonian_elem(det_a, det_a, n_up, n_dn);
  //   pt_energy += pow(sum_a, 2) / (var_energy - H_aa);
  // }
  if (mpi.id == 0) printf("PT energy correction: %.10f\n", pt_energy);
}
  
}
