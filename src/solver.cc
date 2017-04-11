#include "solver.h"

#include <algorithm>
#include <boost/functional/hash.hpp>
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/utility.hpp>
#include <chrono>
#include <clocale>
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
#include "status.h"
#include "types.h"
#include "wavefunction.h"

namespace hci {

Solver::Solver() {
  mpi.id = mpi.world.rank();
  mpi.n = mpi.world.size();
  setlocale(LC_NUMERIC, "");
}

// Main solve procedure.
void Solver::solve() {
  std::string proc_name = boost::mpi::environment::processor_name();
  printf("Proc %d running on %s\n", mpi.id, proc_name.c_str());
  mpi.world.barrier();

  if (mpi.id == 0) printf("%s Begin solving.\n", Status::time());
  std::ifstream config_file("CONFIG");
  if (!config_file) {
    if (mpi.id == 0) printf("CONFIG file not found.\n");
    return;
  }
  read_config(config_file);
  setup();
  load_wavefunction(wave_filename);
  pt_det(eps_pt);
  const double cor_energy = var_energy + pt_energy - hf_energy;
  if (mpi.id == 0) printf("Correlation energy: %.10f\n", cor_energy);
}

void Solver::load_wavefunction(const std::string& filename) {
  int n_in;
  double coef;
  const Det zero_det(n_orbs);
  int orb;

  // Read header line.
  std::ifstream wave_file(filename.c_str());
  if (!wave_file) {
    if (mpi.id == 0) printf("Wave file %s not found.\n", filename.c_str());
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
    printf("%s Loaded var dets (%'d)\n", Status::time(), wf.size());
  }
}

// Deterministic 2nd-order purterbation.
void Solver::pt_det(const double eps_pt) {
  mpi.world.barrier();
  auto begin = std::chrono::high_resolution_clock::now();
  if (mpi.id == 0) {
    printf("%s Start PT w/ %d procs.\n", Status::time(), mpi.n);
  }

  // Save variational dets into hash set.
  std::unordered_set<Det, boost::hash<Det>> var_dets_set;
  const auto& var_dets = wf.get_dets();
  for (const auto& det : var_dets) {
    var_dets_set.insert(det);
  }
  var_dets_set.rehash(var_dets_set.size() * 10);

  // Determine number of hash buckets.
  // if (mpi.id == 0) Status::print("Estimating number of PT dets.");
  const int n = wf.size();
  unsigned long long hash_buckets_local = 0, hash_buckets = 0;
  const auto& var_coefs = wf.get_coefs();
  auto it_det = var_dets.begin();
  auto it_coef = var_coefs.begin();
  const int sample_interval = std::max(n / 500, 100);
  for (int i = 0; i < n; i++) {
    const auto& det_i = *it_det++;
    const double coef_i = *it_coef++;
    if ((i % (sample_interval * mpi.n)) != sample_interval * mpi.id) continue;
    const auto& connected_dets =
        find_connected_dets(det_i, eps_pt / fabs(coef_i));
    hash_buckets_local += connected_dets.size();
  }
  hash_buckets_local *= sample_interval;
  all_reduce(
      mpi.world,
      hash_buckets_local,
      hash_buckets,
      std::plus<unsigned long long>());
  if (mpi.id == 0) printf("Estimated PT dets: %'llu\n", hash_buckets);
  hash_buckets *= 2;

  // Accumulate pt_sum for each det_a.
  std::pair<std::pair<EncodeType, EncodeType>, double> skeleton;
  skeleton.first = var_dets.front().encode();
  BigUnorderedMap<
      std::pair<EncodeType, EncodeType>,
      double,
      boost::hash<std::pair<EncodeType, EncodeType>>>
      pt_sums(mpi.world, skeleton, 200);
  pt_sums.reserve(hash_buckets);
  hash_buckets = pt_sums.bucket_count();
  if (mpi.id == 0) {
    printf(
        "%s Reserved %'llu total hash buckets.\n",
        Status::time(),
        hash_buckets);
  }
  it_det = var_dets.begin();
  it_coef = var_coefs.begin();
  std::vector<Det> var_dets_shrinked;
  std::vector<double> var_coefs_shrinked;
  var_dets_shrinked.reserve(n / mpi.n + 1);
  var_coefs_shrinked.reserve(n / mpi.n + 1);
  for (int i = 0; i < n; i++) {
    const auto& det_i = *it_det++;
    const double coef_i = *it_coef++;
    if (i % mpi.n != mpi.id) continue;
    var_dets_shrinked.push_back(det_i);
    var_coefs_shrinked.push_back(coef_i);
  }
  wf.clear();
  int progress = 10;
  const int local_n = static_cast<int>(var_dets_shrinked.size());
  for (int i = 0; i < local_n; i++) {
    const auto& det_i = var_dets_shrinked[i];
    const double coef_i = var_coefs_shrinked[i];
    const auto& connected_dets =
        find_connected_dets(det_i, eps_pt / fabs(coef_i));
    for (const auto& det_a : connected_dets) {
      if (var_dets_set.count(det_a) == 1) continue;
      const double H_ai = get_hamiltonian_elem(det_i, det_a, n_up, n_dn);
      if (fabs(H_ai) < Constants::EPSILON) continue;
      const double term = H_ai * coef_i;
      pt_sums.async_inc(det_a.encode(), term);
    }
    if ((i + 1) * 100 >= local_n * progress && mpi.id == 0) {
      const auto& local_map = pt_sums.get_local_map();
      const std::size_t local_size = local_map.size();
      const double load_factor = local_map.load_factor();
      printf("%s ", Status::time());
      printf("MASTER: Progress: %d%% (%d/%d) ", progress, i, local_n);
      printf("Local PT dets: %'lu, hash load: %.2f\n", local_size, load_factor);
      progress += 10;
    }
  }
  pt_sums.complete_async_incs();

  unsigned long long n_pt_dets = pt_sums.size();
  if (mpi.id == 0) {
    printf(
        "%s Accumalation done, total PT dets: %'llu\n",
        Status::time(),
        n_pt_dets);
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

  if (mpi.id == 0) {
    using namespace std::chrono;
    printf("%s Sum done. PT correction: %.10f\n", Status::time(), pt_energy);
    auto end = high_resolution_clock::now();
    auto pt_time = duration_cast<duration<double>>(end - begin).count();
    printf("%s Total PT time: %.3f\n", Status::time(), pt_time);
  }
}
}
