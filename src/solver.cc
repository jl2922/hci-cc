#include "solver.h"

#include <algorithm>
#include <boost/functional/hash.hpp>
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/utility.hpp>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <list>
#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include "big_unordered_map.h"
#include "constants.h"
#include "det.h"
#include "types.h"
#include "wavefunction.h"

namespace hci {

Solver::Solver() {
  mpi.id = mpi.world.rank();
  mpi.n = mpi.world.size();
}

// Main solve procedure.
void Solver::solve() {
  std::ifstream config_file("CONFIG");
  if (!config_file) {
    if (mpi.id == 0) std::cout << "CONFIG file not found." << std::endl;
    return;
  }
  read_config(config_file);
  setup();
  load_wavefunction(wave_filename);
  pt_det(eps_pt);
  if (mpi.id == 0) {
    printf("Correlation energy: %.10f\n", var_energy + pt_energy - hf_energy);
  }
}

void Solver::load_wavefunction(const std::string& filename) {
  int n_in;
  double coef;
  const Det zero_det(n_orbs);
  int orb;

  // Read header line.
  std::ifstream wave_file(filename.c_str());
  if (!wave_file) {
    if (mpi.id == 0) {
      std::cout << "Wave file (" << filename << ") not found." << std::endl;
    }
    exit(0);
  }
  wave_file >> n_in >> hf_energy >> var_energy;

  // Read each coef and det.
  for (int i = 0; i < n_in; i++) {
    wave_file >> coef;
    Det& det = wf.append_det(zero_det, coef);
    for (int j = 0; j < n_up; j++) {
      wave_file >> orb;
      det.up.set_orb(orb - 1, true);
    }
    for (int j = 0; j < n_dn; j++) {
      wave_file >> orb;
      det.dn.set_orb(orb - 1, true);
    }
  }

  mpi.world.barrier();
  if (mpi.id == 0) {
    printf("Loaded wavefunction w/ %d dets.\n", wf.size());
    fflush(stdout);
  }
}

// Deterministic 2nd-order purterbation.
void Solver::pt_det(const double eps_pt) {
  if (mpi.id == 0) {
    printf("Performing PT with %d procs.\n", mpi.n);
    fflush(stdout);
  }
  auto begin = std::chrono::high_resolution_clock::now();

  // Save variational dets into hash set.
  std::unordered_set<Det, boost::hash<Det>> var_dets_set;
  const auto& var_dets = wf.get_dets();
  for (const auto& det : var_dets) {
    var_dets_set.insert(det);
  }
  var_dets_set.rehash(var_dets_set.size() * 10);

  // Determine number of hash buckets.
  const int n = wf.size();
  int hash_buckets_local = 0, hash_buckets = 0;
  const auto& var_coefs = wf.get_coefs();
  auto it_det = var_dets.begin();
  auto it_coef = var_coefs.begin();
  const int sample_interval = std::min(n / 500, 100);
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
  std::pair<std::pair<EncodeType, EncodeType>, double> skeleton;
  skeleton.first = var_dets.front().encode();
  BigUnorderedMap<
      std::pair<EncodeType, EncodeType>,
      double,
      boost::hash<std::pair<EncodeType, EncodeType>>>
      pt_sums(mpi.world, skeleton, 200);
  pt_sums.reserve(hash_buckets);
  if (mpi.id == 0) {
    printf("Reserve hash map with %d buckets.\n", hash_buckets);
    fflush(stdout);
  }
  it_det = var_dets.begin();
  it_coef = var_coefs.begin();
  int progress = 10;
  for (int i = 0; i < n; i++) {
    const auto& det_i = *it_det++;
    const double coef_i = *it_coef++;
    if (i % mpi.n != mpi.id) continue;
    const auto& connected_dets =
        find_connected_dets(det_i, eps_pt / fabs(coef_i));
    for (const auto& det_a : connected_dets) {
      if (var_dets_set.count(det_a) == 1) continue;
      const double H_ai = get_hamiltonian_elem(det_i, det_a, n_up, n_dn);
      if (fabs(H_ai) < Constants::EPSILON) continue;
      const double term = H_ai * coef_i;
      pt_sums.async_inc(det_a.encode(), term);
    }
    if ((i + 1) * 100 >= n * progress && mpi.id == 0) {
      const auto& local_map = pt_sums.get_local_map();
      printf(
          "MASTER: Progress: %d%% (%d/%d), PT dets: %lu, hash load: %.2f\n",
          progress,
          i,
          n,
          local_map.size(),
          local_map.load_factor());
      fflush(stdout);
      progress += 10;
    }
  }
  pt_sums.complete_async_incs();

  std::size_t n_pt_dets = pt_sums.size();
  if (mpi.id == 0) {
    printf("Number of PT dets: %lu\n", n_pt_dets);
    fflush(stdout);
  }

  // Accumulate contribution from each det_a to the pt_energy.
  pt_energy = 0.0;
  Det det_a;
  double pt_energy_local = 0.0;
  for (const auto& kv : pt_sums.get_local_map()) {
    det_a.decode(kv.first, n_orbs);
    const double sum_a = kv.second;
    const double H_aa = get_hamiltonian_elem(det_a, det_a, n_up, n_dn);
    pt_energy_local += pow(sum_a, 2) / (var_energy - H_aa);
  }
  reduce(mpi.world, pt_energy_local, pt_energy, std::plus<double>(), 0);
  if (mpi.id == 0) printf("PT energy correction: %.10f\n", pt_energy);
  auto end = std::chrono::high_resolution_clock::now();
  if (mpi.id == 0) {
    std::chrono::duration<double> pt_time = end - begin;
    std::cout << "PT Time: " << pt_time.count() << " s" << std::endl;
  }
}
}
