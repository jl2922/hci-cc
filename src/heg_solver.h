#ifndef HCI_HEG_SOLVER_H_
#define HCI_HEG_SOLVER_H_

#include <array>
#include <valarray>
#include <vector>
#include <unordered_map>
#include <utility>
#include <boost/functional/hash.hpp>
#include "array_math.h"
#include "det.h"
#include "solver.h"
#include "types.h"
#include "wavefunction.h"

namespace hci {

struct HEGData {
  double cell_length;
  double k_unit;
  double r_cutoff;
  double r_s;
  std::vector<Int3> k_vectors;
  std::unordered_map<Int3, int, boost::hash<Int3>> orb_lut;
  std::unordered_map<Int3, std::vector<Int3Double>, boost::hash<Int3>> same_spin_queue;
  std::vector<Int3Double> opposite_spin_queue;
  int n_diff;
  int n_diff_offset;
  int n_max;
};

class HEGSolver: public Solver {
  public:
  protected:
    Wavefunction find_connected_dets(const Det&, const double);
    double get_hamiltonian_elem(const Det&, const Det&) const;
    void setup();
    HEGData heg;
  private:
    int find_orb_id(const std::array<int, 3>&);
    int get_gamma_exp(
        const SpinDet&,
        const std::vector<int>&) const;
    void generate_k_vectors();
    void generate_hci_queue();
    void generate_orb_lut();
    double get_abs_hamiltonian_by_pqrs(
        const int, const int, const int, const int) const;
};

}

#endif // HCI_HEG_SOLVER_H_
