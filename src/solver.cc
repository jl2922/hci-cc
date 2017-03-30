#include "solver.h"

#include <unordered_map>
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
  double pt_energy = 0.0;
  std::unordered_map<hci::Det, int, boost::hash<hci::Det>> var_dets;
  int n = this->wf.n;

  for (int i = 0; i < n; i++) {
    auto& det_i = this->wf.get_det(i);
    var_dets[det_i] = i;
  }
  for (int i = 0; i < n; i++) {
    auto& det_i = this->wf.get_det(i);
    auto coef_i = this->wf.get_coef(i);
    const auto& connected_dets =
        this->find_connected_dets(det_i, eps_pt / abs(coef_i));
  }
}
  
}