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
#include "timer.h"
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

  if (mpi.id == 0) printf("%s Begin solving.\n", Timer::str());
  std::ifstream config_file("CONFIG");
  if (!config_file) {
    if (mpi.id == 0) printf("CONFIG file not found.\n");
    return;
  }
  read_config(config_file);
  setup();
  load_wavefunction(wave_filename);
  pt(eps_pt);
  const double cor_energy = var_energy + pt_energy - hf_energy;
  if (mpi.id == 0) {
    printf("\n[Summary]\n");
    printf("HF Energy: %.10f eV\n", hf_energy);
    printf("Variational Energy: %.10f eV\n", var_energy);
    printf("PT Correction: %.10f eV\n", pt_energy);
    printf("Correlation energy: %.10f eV\n", cor_energy);
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
    printf("%s Loaded var dets (%'d)\n", Timer::str(), wf.size());
  }
}

unsigned long long Solver::estimate_n_pt_dets() {
  const int n = wf.size();
  unsigned long long estimation_local = 0, estimation = 0;
  auto it_det = wf.get_dets().begin();
  auto it_coef = wf.get_coefs().begin();
  const int sample_interval = std::max(n / 500, 100);
  for (int i = 0; i < n; i++) {
    const auto& det_i = *it_det++;
    const double coef_i = *it_coef++;
    if ((i % (sample_interval * mpi.n)) != sample_interval * mpi.id) continue;
    const auto& connected_dets =
        find_connected_dets(det_i, eps_pt / fabs(coef_i));
    estimation_local += connected_dets.size();
  }
  estimation_local *= sample_interval;
  all_reduce(
      mpi.world, estimation_local, estimation, std::plus<unsigned long long>());
  return estimation;
}

template <class T>
std::vector<T> Solver::get_local_portion(const std::list<T>& all) {
  std::vector<T> local;
  auto it = all.begin();
  local.reserve(all.size() / mpi.n + 1);
  for (int i = 0; i < static_cast<int>(all.size()); i++) {
    const auto& item = *it++;
    if (i % mpi.n != mpi.id) continue;
    local.push_back(item);
  }
  return local;
}

// Deterministic 2nd-order purterbation.
void Solver::pt(const double eps_pt) {
  mpi.world.barrier();
  auto begin = std::chrono::high_resolution_clock::now();
  if (mpi.id == 0) printf("%s Start PT w/ %d procs.\n", Timer::str(), mpi.n);

  // Save variational dets into hash set.
  std::unordered_set<Det, boost::hash<Det>> var_dets_set;
  for (const auto& det : wf.get_dets()) var_dets_set.insert(det);
  var_dets_set.rehash(var_dets_set.size() * 10);

  // Determine number of hash buckets.
  unsigned long long n_pt_dets_estimate = estimate_n_pt_dets();
  if (mpi.id == 0) printf("Estimated PT dets: %'llu\n", n_pt_dets_estimate);

  // Setup hash table.
  std::pair<std::pair<EncodeType, EncodeType>, double> skeleton;
  skeleton.first = wf.get_dets().front().encode();
  BigUnorderedMap<
      std::pair<EncodeType, EncodeType>,
      double,
      boost::hash<std::pair<EncodeType, EncodeType>>>
      pt_sums(mpi.world, skeleton);
  pt_sums.reserve(n_pt_dets_estimate * 2);
  unsigned long long hash_buckets = pt_sums.bucket_count();
  if (mpi.id == 0) {
    printf("%s ", Timer::str());
    printf("Reserved %'llu total hash buckets.\n", hash_buckets);
  }

  // Shrink variational dets.
  const auto& var_dets_shrinked = get_local_portion(wf.get_dets());
  const auto& var_coefs_shrinked = get_local_portion(wf.get_coefs());
  wf.clear();

  // Accumulate sums.
  int progress = 10;
  std::size_t local_n = var_dets_shrinked.size();
  for (std::size_t i = 0; i < local_n; i++) {
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
      printf("%s ", Timer::str());
      printf("MASTER: Progress: %d%% (%lu/%lu) ", progress, i, local_n);
      printf("Local PT dets: %'lu, hash load: %.2f\n", local_size, load_factor);
      progress += 10;
    }
  }
  pt_sums.complete_async_incs();

  unsigned long long n_pt_dets = pt_sums.size();
  if (mpi.id == 0) printf("%s Total PT dets: %'llu\n", Timer::str(), n_pt_dets);

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
    printf("%s PT correction: %.10f eV\n", Timer::str(), pt_energy);
    auto end = high_resolution_clock::now();
    auto pt_time = duration_cast<duration<double>>(end - begin).count();
    printf("%s Total PT time: %.3f s\n", Timer::str(), pt_time);
  }
}
}
