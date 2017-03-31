#include "solver.h"

#include <cstdio>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include "det.h"
#include "wavefunction.h"

namespace hci {

void Solver::solve() {
  this->setup();
  this->wf.load("test/WAVE_small");
  this->pt_det(1.0e-5);
}

void Solver::pt_det(const double eps_pt) {
  std::unordered_set<hci::Det, boost::hash<hci::Det>> var_dets;
  std::list<hci::Det> pt_dets;
  std::unordered_map<hci::Det, double, boost::hash<hci::Det>> pt_sums;
  int n = this->wf.n;

  for (int i = 0; i < n; i++) {
    auto& det_i = this->wf.get_det(i);
    var_dets.insert(det_i);
  }
  for (int i = 0; i < n; i++) {
    const auto& det_i = this->wf.get_det(i);
    const double coef_i = this->wf.get_coef(i);
    const auto& connected_dets =
        this->find_connected_dets(det_i, eps_pt / abs(coef_i));
    for (int a = 0; a < connected_dets.n; a++) {
      const auto& det_a = connected_dets.get_det(a);
      if (var_dets.count(det_a) == 1) continue;
      const double H_ai = this->get_hamiltonian_elem(det_i, det_a);
      if (abs(H_ai) < 1.0e-10) continue;
      const double term = H_ai * coef_i;
      if (pt_sums.count(det_a) == 1) {
        pt_sums[det_a] += term; 
      } else {
        pt_dets.push_back(det_a);
        pt_sums[det_a] = term;
      }
    }
  }

  double pt_energy = 0.0;
  for (auto it = pt_dets.begin(); it != pt_dets.end(); it++) {
    auto& det_a = *it;
    const double H_aa = this->get_hamiltonian_elem(det_a, det_a);
    pt_energy += pow(pt_sums[det_a], 2) / (this->var_energy - H_aa);
  }
  this->pt_energy = pt_energy;
}
  
}