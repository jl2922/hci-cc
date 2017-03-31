#include "heg_solver.h"

#include <cmath>
#include "det.h"
#include "solver.h"
#include "wavefunction.h"

namespace hci {

double HEGSolver::get_hamiltonian_elem(const Det& det_pq, const Det& det_rs) const {
  auto& heg = this->heg;
  auto& k_vectors = heg.k_vectors;
  double k_unit = heg.k_unit;
  double one_over_pi_l = 1.0 / (M_PI * heg.cell_length);
  double H = 0.0;

  if (det_pq == det_rs) {
    // Diagonal elements.
    int n_occ_up = det_pq.up.get_n_elec();
    int n_occ_dn = det_pq.dn.get_n_elec();
    const auto& occ_pq_up = det_pq.up.get_elec_orbs(n_occ_up);
    const auto& occ_pq_dn = det_pq.dn.get_elec_orbs(n_occ_dn);

    // One electron operator.
    for (int i = 0; i < n_occ_up; i++) {
      int p = occ_pq_up[i];
      double sum = 0.0;
      for (int k = 0; k < 3; k++) {
        sum += pow(k_vectors[p][k] * k_unit, 2) * 0.5;
      }
      H += sum;
    }
    for (int i = 0; i < n_occ_dn; i++) {
      int p = occ_pq_dn[i];
      double sum = 0.0;
      for (int k = 0; k < 3; k++) {
        sum += pow(k_vectors[p][k] * k_unit, 2) * 0.5;
      }
      H += sum;
    }

    // Two electrons operator.
    for (int i = 0; i < n_occ_up; i++) {
      int p = occ_pq_up[i];
      for (int j = i + 1; j < n_occ_up; j++) {
        int q = occ_pq_up[j];
        double sum = 0.0;
        for (int k = 0; k < 3; k++) {
          sum += pow(k_vectors[p][k] - k_vectors[q][k], 2);
        }
        H -= one_over_pi_l / sum;
      }
    }
    for (int i = 0; i < n_occ_dn; i++) {
      int p = occ_pq_dn[i];
      for (int j = i + 1; j < n_occ_dn; j++) {
        int q = occ_pq_dn[j];
        double sum = 0.0;
        for (int k = 0; k < 3; k++) {
          sum += pow(k_vectors[p][k] - k_vectors[q][k], 2);
        }
        H -= one_over_pi_l / sum;
      }
    }
  } else {
    // Off-diagonal elements.
    static Det det_eor;
    det_eor.from_eor(det_pq, det_rs);
    int n_eor_up = det_eor.up.get_n_elec();
    int n_eor_dn = det_eor.dn.get_n_elec();
    if (n_eor_up + n_eor_dn != 4) return 0.0;
    const auto& eor_up_set_bits = det_eor.up.get_elec_orbs(n_eor_up);
    const auto& eor_dn_set_bits = det_eor.dn.get_elec_orbs(n_eor_dn);
    int k_change[3] = {0, 0, 0};
    bool k_p_set = false, k_q_set = false, k_s_set = false;
    int orb_p = 0, orb_q = 0, orb_s = 0;

    for (int i = 0; i < n_eor_up; i++) {
      int orb_i = eor_up_set_bits[i];
      if (det_pq.up.get_orb(orb_i)) {
        for (int k = 0; k < 3; k++) k_change[k] -= k_vectors[orb_i][k];
        if (!k_p_set) {
          orb_p = orb_i;
          k_p_set = true;
        }
      } else {
        for (int k = 0; k < 3; k++) k_change[k] += k_vectors[orb_i][k];
        if (!k_q_set) {
          orb_q = orb_i;
          k_q_set = true;
        } else {
          orb_s = orb_i;
          k_s_set = true;
        }
      }
    }
    for (int i = 0; i < n_eor_dn; i++) {
      int orb_i = eor_dn_set_bits[i];
      if (det_pq.dn.get_orb(orb_i)) {
        for (int k = 0; k < 3; k++) k_change[k] -= k_vectors[orb_i][k];
        if (!k_p_set) {
          orb_p = orb_i;
          k_p_set = true;
        }
      } else {
        for (int k = 0; k < 3; k++) k_change[k] += k_vectors[orb_i][k];
        if (!k_q_set) {
          orb_q = orb_i;
          k_q_set = true;
        } else {
          orb_s = orb_i;
          k_s_set = true;
        }
      }
    }
    
    for (int k = 0; k < 3; k++) {
      if (k_change[k] != 0) return 0.0; // Momentum not conserved.
    }

    
  }

  return H;
}

void HEGSolver::setup() {

}

Wavefunction HEGSolver::find_connected_dets(const Det& det, const double eps) const {
  Wavefunction wf;
  int n_pq_pairs = 0;
  int n_rs_pairs = 0;

  if (this->max_abs_H < eps) {
    return wf;
  }

  for (int i = 0; i < n_pq_pairs; i++) {
    for (int j = 0; j < n_rs_pairs; j++) {

    }
  }
  return wf;
}
  
}
