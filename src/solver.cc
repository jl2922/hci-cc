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
#include "constants.h"
#include "det.h"
#include "parallel.h"
#include "types.h"
#include "wavefunction.h"
#include "gps.h"

std::mutex m_screen;

namespace hci {

Solver::Solver() {
  mpi.id = mpi.world.rank();
  mpi.n = mpi.world.size();
  mpi.is_master = (mpi.id == 0);
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
  if (mpi.is_master) printf("Loaded wavefunction w/ %d dets.\n", wf.size());
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
  // Generate variational dets hash set.
  const int n = wf.size();
  const auto& var_dets = wf.get_dets();
  std::unordered_set<Det, boost::hash<Det>> var_dets_set(n * 5);
  for (const auto& det: var_dets) var_dets_set.insert(det);

  // Determine number of hash buckets.
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
  hash_buckets_local *= 2 * sample_interval / mpi.n;
  all_reduce(mpi.world, hash_buckets_local, hash_buckets, std::plus<int>());
  std::unordered_map<Det, double, boost::hash<Det>> pt_sums(hash_buckets);
  if (mpi.id == 0) printf("Hash buckets PP: %lu\n", pt_sums.bucket_count());

  // Accumulate pt_sum for each det_a.
  it_det = var_dets.begin();
  it_coef = var_coefs.begin();
  int progress = 10;
  std::pair<Det, double> proto;
  proto.first = Det(n_orbs);
  proto.second = 0.0;
  const int BUF_SIZE = 100;
  std::array<std::pair<Det, double>, BUF_SIZE> buf_send, buf_recv;
  std::array<boost::mpi::content, BUF_SIZE> contents_send, contents_recv;
  for (int i = 0; i < BUF_SIZE; i++) {
    buf_send[i] = proto;
    boost::mpi::broadcast(mpi.world, boost::mpi::skeleton(buf_send[i]), 0);
    contents_send[i] = boost::mpi::get_content(buf_send[i]);
    buf_recv[i] = proto;
    boost::mpi::broadcast(mpi.world, boost::mpi::skeleton(buf_recv[i]), 0);
    contents_recv[i] = boost::mpi::get_content(buf_recv[i]);
  }
  mpi.world.barrier();
  int n_send = 0, n_recv = 0;
  std::list<boost::mpi::request> reqs;
  for (int i = 0; i < 20; i++) {
    mpi.world.barrier();
    if (mpi.id == 0) printf("i = %d\n", i);
    if (i % mpi.n != mpi.id) {
      
    } else {
      const auto& det_i = *it_det++;
      const double coef_i = *it_coef++;
      const auto& connected_dets =
          find_connected_dets(det_i, eps_pt / fabs(coef_i));
      int buf_id = 0;
      for (const auto& det_a: connected_dets) {
        if (var_dets_set.count(det_a) == 1) continue;
        const double H_ai = get_hamiltonian_elem(det_i, det_a, n_up, n_dn);
        if (fabs(H_ai) < Constants::EPSILON) continue;
        const double term = H_ai * coef_i;
        int target = static_cast<int>(abs(hash_value(det_a))) % mpi.n;
        if (target == mpi.id) {
          pt_sums[det_a] += term;
        } else {
          buf_send[buf_id].first = det_a;
          buf_send[buf_id].second = term;
          reqs.push_front(mpi.world.isend(target, TAG_PT_TERM, contents_send[buf_id]));
          n_send++;
          buf_id++;
          if (buf_id == BUF_SIZE) {
            wait_all(reqs.begin(), reqs.end());
            reqs.clear();
            buf_id = 0;
          }
        }
      }
    }
    printf("n_send: %d\n", n_send);
    wait_all(reqs.begin(), reqs.end());
    reqs.clear();
    int buf_id = 0;
    while (mpi.world.iprobe(boost::mpi::any_source, TAG_PT_TERM)) {
      reqs.push_front(mpi.world.irecv(boost::mpi::any_source, TAG_PT_TERM, contents_recv[buf_id]));
      buf_id++;
      if (buf_id == BUF_SIZE) {
        wait_all(reqs.begin(), reqs.end());
        auto it_req = reqs.begin();
        for (int a = 0; a < BUF_SIZE; a++) {
          it_req++;
          const auto& det_a = buf_recv[a].first;
          n_recv++;
          if (det_a.is_zero()) continue;
          const double term = buf_recv[a].second;
          int target = static_cast<int>(abs(hash_value(det_a))) % mpi.n;
          if (target != mpi.id) {
            std::stringstream ss;
            ss << "From: " << (*it_req).wait().source() << std::endl;
            ss << det_a.up << std::endl;
            ss << det_a.dn << std::endl;
            ss << term << std::endl;
            printf("Wrong delivery. Should go to %d, received by %d.\n%s", target, mpi.id, ss.str().c_str());
            fflush(stdout);
            exit(1);
          }
          pt_sums[det_a] += term;
        }
        reqs.clear();
        buf_id = 0;
      }
    }
    printf("n_recv: %d\n", n_recv);
    wait_all(reqs.begin(), reqs.end());
    // for (const auto& term: buffers) {
    //   const auto& det_a = term.first;
    //   pt_sums[det_a] = term.second;
    // }
    if (reqs.size() != buf_id) {
      printf("Reqs size: %lu, buf id: %d\n", reqs.size(), buf_id);
      fflush(stdout);
      exit(1);
    }
    auto it_req = reqs.begin();
    for (int a = 0; a < buf_id; a++) {
      const auto& det_a = buf_recv[a].first;
      it_req++;
      n_recv++;
      // if (det_a.is_zero()) {
      //   continue;
      // }
      const double term = buf_recv[a].second;
      int target = static_cast<int>(abs(hash_value(det_a))) % mpi.n;
      if (target != mpi.id || det_a.is_zero()) {
        std::stringstream ss;
        ss << "From: " << (*it_req).wait().source() << std::endl;
        ss << det_a.up << std::endl;
        ss << det_a.dn << std::endl;
        ss << term << std::endl;
        printf("Wrong delivery #1. Should go to %d, received by %d.\n%s", target, mpi.id, ss.str().c_str());
        fflush(stdout);
        // exit(1);
      }
      pt_sums[det_a] += term;
    }
    reqs.clear();
    printf("n_recv: %d\n", n_recv);
    // buffers.clear();
    // for (int i = 0; i < mpi.n; i++) {
    //   if (i == mpi.id) continue;
    //   mpi.world.isend(i, TAG_PT_BARRIER, true);
    // }
    // int finished_count = 0;
    // while (finished_count < mpi.n) {
    //   const auto& status = mpi.world.probe();
    //   const int source = status.source();
    //   const auto tag = status.tag();
    //   if (tag == TAG_PT_TERM) {
    //     mpi.world.recv(source, TAG_PT_TERM, skeleton);
    //     printf("Proc %d received a det from proc %d.\n", mpi.id, source);
    //     const auto& det_a = skeleton.first;
    //     pt_sums[det_a] = skeleton.second;
    //   } else if (tag == TAG_PT_BARRIER) {
    //     mpi.world.recv(source, TAG_PT_TERM);
    //     finished_count++;
    //   }
    // }

    if (mpi.is_master && i * 100 >= (n - 1) * progress) {
      printf("Master Progress: %d%% (%d/%d) PT dets: %lu, hash load: %.2f\n",
          progress, i, n, pt_sums.size(), pt_sums.load_factor());
      fflush(stdout);
      progress += 10;
    }
  }
  printf("Process %d finished. # PT dets: %lu\n", mpi.id, pt_sums.size());
  mpi.world.barrier();
  printf("n_send, n_recv: %d, %d\n", n_send, n_recv);
  int buf_id = 0;
  while (mpi.world.iprobe(boost::mpi::any_source, TAG_PT_TERM)) {
    // printf("processe %d received a det.\n", mpi.id);
    // fflush(stdout);
    // buffers.push_back(skeleton);
    reqs.push_front(mpi.world.irecv(boost::mpi::any_source, TAG_PT_TERM, contents_recv[buf_id]));
    buf_id++;
    if (buf_id == BUF_SIZE) {
      wait_all(reqs.begin(), reqs.end());
      reqs.clear();
      for (int a = 0; a < BUF_SIZE; a++) {
        const auto& det_a = buf_recv[a].first;
        n_recv++;
        if (det_a.is_zero()) continue;
        const double term = buf_recv[a].second;
        int target = static_cast<int>(abs(hash_value(det_a))) % mpi.n;
        if (target != mpi.id) {
        //   std::stringstream ss;
        //   ss << det_a.up << std::endl;
        //   ss << det_a.dn << std::endl;
        //   ss << term << std::endl;
          // printf("Wrong delivery. Should go to %d, received by %d.\n%s", target, mpi.id, ss.str().c_str());
          // fflush(stdout);
          exit(1);
        }
        pt_sums[det_a] += term;
      }
      buf_id = 0;
    }
  }
  printf("n_recv: %d\n", n_recv);
  wait_all(reqs.begin(), reqs.end());
  // for (const auto& term: buffers) {
  //   const auto& det_a = term.first;
  //   pt_sums[det_a] = term.second;
  // }
  reqs.clear();
  for (int a = 0; a < buf_id; a++) {
    const auto& det_a = buf_recv[a].first;
    n_recv++;
    if (det_a.is_zero()) continue;
    const double term = buf_recv[a].second;
    int target = static_cast<int>(abs(hash_value(det_a))) % mpi.n;
    if (target != mpi.id) {
      std::stringstream ss;
      ss << det_a.up << std::endl;
      ss << det_a.dn << std::endl;
      ss << term << std::endl;
      // printf("Wrong delivery. Should go to %d, received by %d.\n%s", target, mpi.id, ss.str().c_str());
      fflush(stdout);
      exit(1);
    }
    pt_sums[det_a] += term;
  }
  mpi.world.barrier();
  printf("n_send, n_recv: %d, %d\n", n_send, n_recv);

  int n_pt_dets = 0;
  reduce(mpi.world, static_cast<int>(pt_sums.size()), n_pt_dets, std::plus<int>(), 0);
  if (mpi.is_master) printf("Total # PT dets: %d\n", n_pt_dets);

  // Accumulate contribution from each det_a to the pt_energy.
  double pt_energy_local = 0.0, pt_energy = 0;
  for (auto it = pt_sums.begin(); it != pt_sums.end(); it++) {
    const auto& det_a = it->first;
    const double sum_a = it->second;
    const double H_aa = get_hamiltonian_elem(det_a, det_a, n_up, n_dn);
    pt_energy_local += pow(sum_a, 2) / (var_energy - H_aa);
  }
  reduce(mpi.world, pt_energy_local, pt_energy, std::plus<double>(), 0);
  if (mpi.is_master) printf("Total PT correction: %.10f eV\n", pt_energy);
}
  
}
