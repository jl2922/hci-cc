#include "heg_solver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <forward_list>
#include <fstream>
#include <iostream>
#include <list>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <valarray>
#include <vector>
#include "array_math.h"
#include "big_unordered_map.h"
#include "constants.h"
#include "det.h"
#include "solver.h"
#include "timer.h"
#include "types.h"
#include "wavefunction.h"
namespace hci {

void HEGSolver::read_config(std::ifstream& config_file) {
  config_file >> n_up >> n_dn >> wave_filename >> eps_var >> eps_pt;
  n_elecs = n_up + n_dn;
  config_file >> heg.r_s >> heg.r_cutoff;
  if (mpi.id == 0) {
    printf("n_up: %d, n_dn: %d, n_elecs: %d\n", n_up, n_dn, n_elecs);
    printf("WAVE file: %s\n", wave_filename.c_str());
    printf("eps_var: %.6E\neps_pt: %.6E\n", eps_var, eps_pt);
    printf("r_s: %.6f\n", heg.r_s);
    printf("r_cutoff: %.6f\n", heg.r_cutoff);
  }
}

double HEGSolver::get_abs_hamiltonian_by_pqrs(
    const int p, const int q, const int r, const int s) const {
  if (p == q || r == s || p == r || q == s || p == s || q == r) return 0.0;
  const Det zero_det(n_orbs);
  Det det_pq = zero_det;
  Det det_rs = zero_det;
  det_pq.set_orb(p, n_orbs).set_orb(q, n_orbs);
  det_rs.set_orb(r, n_orbs).set_orb(s, n_orbs);
  return fabs(get_hamiltonian_elem(det_pq, det_rs, -1, -1));
}

// Number of bits set before each changed bit.
int HEGSolver::get_gamma_exp(
    const SpinDet& spin_det,
    const int n_elecs,
    const std::vector<int>& eor) const {
  int gamma_exp = 0;
  int ptr = 0;
  const auto& occ = spin_det.get_elec_orbs(n_elecs);
  for (const int orb_id : eor) {
    if (!spin_det.get_orb(orb_id)) continue;
    while (occ[ptr] < orb_id) ptr++;
    gamma_exp += ptr;
  }
  return gamma_exp;
}

// Calculate hamiltonian between two determinants.
double HEGSolver::get_hamiltonian_elem(
    const Det& det_pq,
    const Det& det_rs,
    const int n_up,
    const int n_dn) const {
  // n_up n_dn are optional inputs for speed up get_elec_orbs.
  auto& k_vectors = heg.k_vectors;
  const double k_unit = heg.k_unit;
  const double one_over_pi_l = 1.0 / (M_PI * heg.cell_length);
  double H = 0.0;

  if (det_pq == det_rs) {
    // Diagonal elements.
    const auto& occ_pq_up = det_pq.up.get_elec_orbs(n_up);
    const auto& occ_pq_dn = det_pq.dn.get_elec_orbs(n_dn);
    const int n_occ_up = occ_pq_up.size();
    const int n_occ_dn = occ_pq_dn.size();

    // One electron operator.
    for (const int p : occ_pq_up) {
      H += sum(square(k_vectors[p] * k_unit) * 0.5);
    }
    for (const int p : occ_pq_dn) {
      H += sum(square(k_vectors[p] * k_unit) * 0.5);
    }

    // Two electrons operator.
    for (int i = 0; i < n_occ_up; i++) {
      const int p = occ_pq_up[i];
      for (int j = i + 1; j < n_occ_up; j++) {
        const int q = occ_pq_up[j];
        H -= one_over_pi_l / sum(square(k_vectors[p] - k_vectors[q]));
      }
    }
    for (int i = 0; i < n_occ_dn; i++) {
      const int p = occ_pq_dn[i];
      for (int j = i + 1; j < n_occ_dn; j++) {
        const int q = occ_pq_dn[j];
        H -= one_over_pi_l / sum(square(k_vectors[p] - k_vectors[q]));
      }
    }
  } else {
    // Off-diagonal elements.
    Det det_eor;
    det_eor.from_eor(det_pq, det_rs);
    const int n_eor_up = det_eor.up.get_n_elecs();
    const int n_eor_dn = det_eor.dn.get_n_elecs();
    if (n_eor_up + n_eor_dn != 4) return 0.0;
    const auto& eor_up_set_bits = det_eor.up.get_elec_orbs(n_eor_up);
    const auto& eor_dn_set_bits = det_eor.dn.get_elec_orbs(n_eor_dn);
    bool k_p_set = false, k_r_set = false;
    int orb_p = 0, orb_r = 0, orb_s = 0;

    // Obtain p, q, s.
    Int3 k_change;
    k_change.fill(0);
    for (const auto& orb_i : eor_up_set_bits) {
      if (det_pq.up.get_orb(orb_i)) {
        k_change -= k_vectors[orb_i];
        if (!k_p_set) {
          orb_p = orb_i;
          k_p_set = true;
        }
      } else {
        k_change += k_vectors[orb_i];
        if (!k_r_set) {
          orb_r = orb_i;
          k_r_set = true;
        } else {
          orb_s = orb_i;
        }
      }
    }
    for (const auto& orb_i : eor_dn_set_bits) {
      if (det_pq.dn.get_orb(orb_i)) {
        k_change -= k_vectors[orb_i];
        if (!k_p_set) {
          orb_p = orb_i;
          k_p_set = true;
        }
      } else {
        k_change += k_vectors[orb_i];
        if (!k_r_set) {
          orb_r = orb_i;
          k_r_set = true;
        } else {
          orb_s = orb_i;
        }
      }
    }

    // Check for momentum conservation.
    if (k_change != 0) return 0.0;

    H = one_over_pi_l / sum(square(k_vectors[orb_p] - k_vectors[orb_r]));
    if (n_eor_up != 2) {
      H -= one_over_pi_l / sum(square(k_vectors[orb_p] - k_vectors[orb_s]));
    }
    int gamma_exp = get_gamma_exp(det_pq.up, n_up, eor_up_set_bits) +
                    get_gamma_exp(det_pq.dn, n_dn, eor_dn_set_bits) +
                    get_gamma_exp(det_rs.up, n_up, eor_up_set_bits) +
                    get_gamma_exp(det_rs.dn, n_dn, eor_dn_set_bits);
    if ((gamma_exp & 1) == 1) H = -H;
  }
  return H;
}

// Setup HEG environment and hci queues.
void HEGSolver::setup() {
  // Basic quantities.
  const double density = 3.0 / (4.0 * M_PI * pow(heg.r_s, 3));
  const double cell_length = pow(n_elecs / density, 1.0 / 3);
  const double k_unit = 2 * M_PI / cell_length;
  heg.cell_length = cell_length;
  heg.k_unit = k_unit;

  generate_k_vectors();
  n_orbs = heg.k_vectors.size();
  if (mpi.id == 0) {
    printf("%s K vectors generated. ", Timer::str());
    printf("Number of spin orbitals: %d\n", n_orbs * 2);
  }

  generate_orb_lut();
  generate_hci_queue();
}

// Find determinants connected to a det passed in.
// Return as a Wavefunction object with coefs all equal zeroes.
std::list<Det> HEGSolver::find_connected_dets(
    const Det& det, const double eps) {
  std::list<Det> connected_dets;
  connected_dets.push_back(det);
  if (max_abs_H < eps) return connected_dets;
  const auto& k_vectors = heg.k_vectors;

  // Get pq pairs.
  std::forward_list<IntPair> pq_pairs;
  const auto& occ_up = det.up.get_elec_orbs(n_up);
  const auto& occ_dn = det.dn.get_elec_orbs(n_dn);
  for (int i = 0; i < n_up; i++) {
    for (int j = i + 1; j < n_up; j++) {
      pq_pairs.push_front(IntPair(occ_up[i], occ_up[j]));
    }
  }
  for (int i = 0; i < n_dn; i++) {
    for (int j = i + 1; j < n_dn; j++) {
      pq_pairs.push_front(IntPair(occ_dn[i] + n_orbs, occ_dn[j] + n_orbs));
    }
  }
  for (int i = 0; i < n_up; i++) {
    for (int j = 0; j < n_dn; j++) {
      pq_pairs.push_front(IntPair(occ_up[i], occ_dn[j] + n_orbs));
    }
  }

  for (const auto& pq_pair : pq_pairs) {
    const int p = pq_pair.first;
    const int q = pq_pair.second;

    // Get rs pairs.
    int pp = p, qq = q;
    if (p >= n_orbs && q >= n_orbs) {
      pp -= n_orbs;
      qq -= n_orbs;
    } else if (p < n_orbs && q >= n_orbs && p > q - n_orbs) {
      pp = q - n_orbs;
      qq = p + n_orbs;
    }
    bool same_spin = false;
    std::vector<Int3Double>* items_ptr;
    if (pp < n_orbs && qq < n_orbs) {
      same_spin = true;
      const auto& diff_pq = k_vectors[qq] - k_vectors[pp];
      items_ptr = &(heg.same_spin_queue[diff_pq]);
    } else {
      items_ptr = &(heg.opposite_spin_queue);
    }
    const auto& items = *items_ptr;
    int qs_offset = 0;
    if (!same_spin) qs_offset = n_orbs;
    for (const auto& item : items) {
      if (item.second < eps) break;
      const auto& diff_pr = item.first;
      int r = find_orb_id(diff_pr + k_vectors[pp]);
      if (r < 0) continue;
      int s =
          find_orb_id(k_vectors[pp] + k_vectors[qq - qs_offset] - k_vectors[r]);
      if (s < 0) continue;
      if (same_spin && s < r) continue;
      s += qs_offset;
      if (p >= n_orbs && q >= n_orbs) {
        r += n_orbs;
        s += n_orbs;
      } else if (p < n_orbs && q >= n_orbs && p > q - n_orbs) {
        const int tmp = s;
        s = r + n_orbs;
        r = tmp - n_orbs;
      }

      // Test whether pqrs is a valid excitation for det.
      if (det.get_orb(r, n_orbs) || det.get_orb(s, n_orbs)) continue;
      connected_dets.push_back(det);
      Det& new_det = connected_dets.back();
      new_det.set_orb(p, n_orbs, false).set_orb(q, n_orbs, false);
      new_det.set_orb(r, n_orbs).set_orb(s, n_orbs);
    }
  }
  return connected_dets;
}

int HEGSolver::find_orb_id(const std::array<int, 3>& k_vector) {
  auto& orb_lut = heg.orb_lut;
  if (orb_lut.count(k_vector) == 1) return orb_lut[k_vector];
  return -1;
}

// Generate k vectors in ascending order of magnitude, and then x, y, z.
void HEGSolver::generate_k_vectors() {
  // Get all valid k vectors and magnitudes.
  const int n_max = floor(heg.r_cutoff);
  heg.n_max = n_max;
  std::vector<std::array<int, 3>> temp_k_vectors;
  std::vector<int> temp_k_length2;
  for (int i = -n_max; i <= n_max; i++) {
    for (int j = -n_max; j <= n_max; j++) {
      for (int k = -n_max; k <= n_max; k++) {
        const int length2 = i * i + j * j + k * k;
        if (length2 > pow(heg.r_cutoff, 2)) continue;
        const Int3 temp_k_vector({i, j, k});
        temp_k_vectors.push_back(temp_k_vector);
        temp_k_length2.push_back(length2);
      }
    }
  }

  // Sort.
  const int n_orbs = temp_k_vectors.size();
  std::vector<int> order;
  for (int i = 0; i < n_orbs; i++) order.push_back(i);
  std::stable_sort(
      order.begin(),
      order.end(),
      [&temp_k_length2](const int& a, const int& b) -> bool {
        return temp_k_length2[a] < temp_k_length2[b];
      });

  // Push to k_vectors.
  for (int i = 0; i < n_orbs; i++) {
    heg.k_vectors.push_back(temp_k_vectors[order[i]]);
  }
}

// Generate lookup table from a k_vector to its index.
void HEGSolver::generate_orb_lut() {
  for (int i = 0; i < n_orbs; i++) {
    heg.orb_lut[heg.k_vectors[i]] = i;
  }
}

// Generate core HCI queue for fast connected determinants finding.
void HEGSolver::generate_hci_queue() {
  max_n_rs_pairs = 0;
  max_abs_H = 0.0;

  generate_same_spin_hci_queue();
  generate_opposite_spin_hci_queue();
}

void HEGSolver::generate_same_spin_hci_queue() {
  auto& k_vectors = heg.k_vectors;

  // Generate k_diffs.
  std::unordered_set<Int3, boost::hash<Int3>> k_diffs_set;
  std::vector<Int3> k_diffs;
  for (int p = 0; p < n_orbs; p++) {
    for (int q = 0; q < n_orbs; q++) {
      if (p == q) continue;
      const auto& diff_pq = k_vectors[q] - k_vectors[p];
      if (k_diffs_set.count(diff_pq) == 1) continue;
      k_diffs.push_back(diff_pq);
      k_diffs_set.insert(diff_pq);
    }
  }
  k_diffs_set.clear();

  // Generate neighbor_diffs;
  std::vector<std::array<double, 3>> neighbor_diffs;
  for (double i = -2; i <= 2; i++) {
    for (double j = -2; j <= 2; j++) {
      for (double k = -2; k <= 2; k++) {
        neighbor_diffs.push_back({i, j, k});
      }
    }
  }
  std::stable_sort(
      neighbor_diffs.begin(),
      neighbor_diffs.end(),
      [](const std::array<double, 3>& a,
         const std::array<double, 3>& b) -> bool { return norm(a) < norm(b); });

  // Process each diff pair.
  int progress = 20;
  int cnt = 0;
  const int total = k_diffs.size();
  for (const auto& diff_pq : k_diffs) {
    for (const auto& diff_pr : k_diffs) {
      const auto& diff_sr =
          diff_pr + diff_pr - diff_pq;  // Momentum conservation.
      if (norm(diff_sr) > heg.r_cutoff * 2 + Constants::EPSILON) continue;
      const auto& p_point_theoretical = diff_pq * (-0.5);
      for (const auto& neighbor_diff : neighbor_diffs) {
        const auto& p_point_guess = round(p_point_theoretical + neighbor_diff);
        const int p = find_orb_id(p_point_guess);
        if (p == -1) continue;
        const int q = find_orb_id(k_vectors[p] + diff_pq);
        if (q < p) continue;
        const int r = find_orb_id(k_vectors[p] + diff_pr);
        if (r == -1) continue;
        const int s = find_orb_id(k_vectors[p] + k_vectors[q] - k_vectors[r]);
        if (s < r) continue;
        const double H_abs = get_abs_hamiltonian_by_pqrs(p, q, r, s);
        if (H_abs < Constants::EPSILON) continue;
        const auto& item = Int3Double(diff_pr, H_abs);
        heg.same_spin_queue[diff_pq].push_back(item);
        break;
      }
    }
    cnt++;
    if (cnt >= progress * total / 100 && mpi.id == 0) {
      printf("%s HCI same spin queue generation progress: ", Timer::str());
      printf("%d %% (%d / %lu)\n", progress, cnt, k_diffs.size());
      progress += 20;
    }
  }

  // Sort.
  for (auto& kv : heg.same_spin_queue) {
    auto& items = kv.second;
    std::sort(
        items.begin(),
        items.end(),
        [](const Int3Double& a, const Int3Double& b) -> bool {
          return a.second > b.second;
        });
    max_abs_H = std::max(max_abs_H, items.front().second);
    const int n_items = static_cast<int>(items.size());
    max_n_rs_pairs = std::max(max_n_rs_pairs, n_items);
  }
}

void HEGSolver::generate_opposite_spin_hci_queue() {
  auto& k_vectors = heg.k_vectors;
  std::unordered_set<Int3, boost::hash<Int3>> opposite_spin_processed;
  for (int p = 0; p < n_orbs; p++) {
    for (int r = 0; r < n_orbs; r++) {
      const auto& diff_pr = k_vectors[r] - k_vectors[p];
      if (opposite_spin_processed.count(diff_pr) == 1) continue;
      opposite_spin_processed.insert(diff_pr);
      for (int q = p; q < n_orbs; q++) {
        const int s = find_orb_id(k_vectors[p] + k_vectors[q] - k_vectors[r]);
        if (s < 0) continue;
        const double H_abs =
            get_abs_hamiltonian_by_pqrs(p, q + n_orbs, r, s + n_orbs);
        if (H_abs < Constants::EPSILON) continue;
        const auto& item = Int3Double(diff_pr, H_abs);
        heg.opposite_spin_queue.push_back(item);
        break;
      }
    }
  }

  // Sort.
  auto& items = heg.opposite_spin_queue;
  std::sort(
      items.begin(),
      items.end(),
      [](const Int3Double& a, const Int3Double& b) -> bool {
        return a.second > b.second;
      });
  max_abs_H = std::max(max_abs_H, items.front().second);
  const int n_items = static_cast<int>(items.size());
  max_n_rs_pairs = std::max(max_n_rs_pairs, n_items);

  if (mpi.id == 0) printf("%s HCI opposite queue generated.\n", Timer::str());
}
}
