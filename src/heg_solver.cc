#include "heg_solver.h"

#include <cmath>
#include <algorithm>
#include <vector>
#include "boost/multi_array.hpp"
#include "det.h"
#include "solver.h"
#include "wavefunction.h"

#define R_S 0.5
#define R_CUTOFF 1.5
#define N_ELEC 14
#define N_UP 7
#define N_DN 7

namespace hci {

double HEGSolver::get_hamiltonian_elem(const Det& det_pq, const Det& det_rs) {
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

    return H;
  } else { // det_pq != det_rs
    // Off-diagonal elements.
    static Det det_eor;
    det_eor.from_eor(det_pq, det_rs);
    int n_eor_up = det_eor.up.get_n_elec();
    int n_eor_dn = det_eor.dn.get_n_elec();
    if (n_eor_up + n_eor_dn != 4) return 0.0;
    const auto& eor_up_set_bits = det_eor.up.get_elec_orbs(n_eor_up);
    const auto& eor_dn_set_bits = det_eor.dn.get_elec_orbs(n_eor_dn);
    int k_change[3] = {0, 0, 0};

    // Obtain pqrs.
    bool k_p_set = false, k_q_set = false;
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
        }
      }
    }

    for (int k = 0; k < 3; k++) {
      if (k_change[k] != 0) return 0.0; // Momentum not conserved.
    }
    double sum = 0.0;
    for (int k = 0; k < 3; k++) {
      sum += pow(k_vectors[orb_p][k] - k_vectors[orb_q][k], 2);
    }
    H = one_over_pi_l / sum;
    if (n_eor_up != 2) {
      double sum = 0.0;
      for (int k = 0; k < 3; k++) {
        sum += pow(k_vectors[orb_p][k] - k_vectors[orb_s][k], 2);
      }
      H -= one_over_pi_l / sum;
    }
    const auto& occ_pq_up = det_pq.up.get_elec_orbs();
    const auto& occ_pq_dn = det_pq.dn.get_elec_orbs();
    const auto& occ_rs_up = det_rs.up.get_elec_orbs();
    const auto& occ_rs_dn = det_rs.dn.get_elec_orbs();
    int gamma_exp =
        get_gamma_exp(det_pq.up, occ_pq_up, eor_up_set_bits, n_eor_up) +
        get_gamma_exp(det_pq.dn, occ_pq_dn, eor_dn_set_bits, n_eor_dn) +
        get_gamma_exp(det_rs.up, occ_rs_up, eor_up_set_bits, n_eor_up) +
        get_gamma_exp(det_rs.dn, occ_rs_dn, eor_dn_set_bits, n_eor_dn);
    if ((gamma_exp & 1) == 1) H = -H;
    return H;
  } // det_pq == det_rs
}

int HEGSolver::get_gamma_exp(
    const SpinDet& spin_det,
    const int* occ,
    const int* eor_set_bits,
    const int n_eor_set_bits) {
  int res = 0;
  int pos = 0;
  for (int i = 0; i < n_eor_set_bits; i++) {
    int orb_i = eor_set_bits[i];
    if (!spin_det.get_orb(orb_i)) continue;
    while (occ[pos + 1] < orb_i) {
      pos++;
    }
    res += pos;
  }
  return res;
}

void HEGSolver::setup() {
  this->heg.r_s = R_S;
  this->heg.r_cutoff = R_CUTOFF;
  this->n_elec = N_ELEC;
  this->n_up = N_UP;
  this->n_dn = N_DN;

  double density = 3.0 / (4.0 * M_PI * pow(this->heg.r_s, 3));
  double cell_length = pow(this->n_elec / density, 1.0 / 3);
  double k_unit = 2 * M_PI / cell_length;
  printf("Cell length: %.10f\n", cell_length);
  printf("K unit: %.10f\n", k_unit);
  this->heg.cell_length = cell_length;
  this->heg.k_unit = k_unit;
  int n_orb = this->generate_k_vectors();
  printf("Number of spin orbitals: %d\n", n_orb * 2);
  this->n_orb = n_orb;

  // call this%generate_orbital_lut()
  // call this%generate_hci_queue()
}

int HEGSolver::generate_k_vectors() {
  int n_max = floor(this->heg.r_cutoff);
  int temp_length = round(pow(2 * n_max + 1, 3));
  this->heg.n_max = n_max;
  auto temp_k_vectors = new int[temp_length][3];
  auto temp_k_length2 = new int[temp_length];
  int n_orbs = 0;
  for (int i = -n_max; i <= n_max; i++) {
    for (int j = -n_max; j <= n_max; j++) {
      for (int k = -n_max; k <= n_max; k++) {
        int length2 = pow(i, 2) + pow(j, 2) + pow(k, 2);
        if (length2 > pow(this->heg.r_cutoff, 2)) continue;
        n_orbs++;
        temp_k_vectors[n_orbs][0] = i;
        temp_k_vectors[n_orbs][1] = j;
        temp_k_vectors[n_orbs][2] = k;
        temp_k_length2[n_orbs] = length2;
      }
    }
  }
  std::vector<int> order(n_orbs);
  for (int i = 0; i < n_orbs; i++) order[i] = i;
  std::stable_sort(
      order.begin(), order.end(),
      [&temp_k_length2](const int& a, const int& b) -> bool {
    return temp_k_length2[a] < temp_k_length2[b];
  });
  printf("K Points:\n");
  for (int i = 0; i < n_orbs; i++) {
    for (int k = 0; k < 3; k++) {
      printf("%.10f ", k_vectors[i][k]);
    }
    printf("\n");
  }
  auto& k_vectors = this->heg.k_vectors;
  k_vectors.reserve(n_orbs);
  do i = 1, n_orbs
    k_vectors.push_back();
    this%heg%k_vectors(:, i) = temp_k_vectors(:, order(i))
    write (6, '(A, I0, 3F15.10)') &
        & 'K-point #', i, this%heg%k_vectors(:, i) * this%heg%k_unit
  end do
  exit(0);
}

Wavefunction HEGSolver::find_connected_dets(const Det& det, const double eps) {
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
