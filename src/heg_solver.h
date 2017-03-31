#ifndef HCI_HEG_SOLVER_H_
#define HCI_HEG_SOLVER_H_

#include <array>
#include <unordered_map>
#include <utility>
#include <boost/functional/hash.hpp>
#include "det.h"
#include "solver.h"
#include "wavefunction.h"

namespace hci {

typedef std::array<int, 3> int_3;
typedef std::array<int, 3> int_4;
typedef std::pair<std::array<int, 3>, double> int_3_double;

struct HEGData {
  double cell_length;
  double k_unit;
  double r_cutoff;
  double r_s;
  std::vector<int_3> k_vectors;
  std::unordered_map<int_3, int, boost::hash<int_3>> orb_lut;
  std::unordered_map<int_4, int_3_double, boost::hash<int_4>> same_spin;
  std::unordered_map<int, int_3_double> opposite_spin;
  int n_diff;
  int n_diff_offset;
  int n_max;
};

class HEGSolver: public Solver {
  public:
  protected:
    Wavefunction find_connected_dets(const Det&, const double) const;
    double get_hamiltonian_elem(const Det&, const Det&) const;
    void setup();
    HEGData heg;
  private:
    int find_orb_id(const std::array<int, 3>&);
    int get_gamma_exp(const SpinDet&, const int*, const int*, const int);
    int generate_k_vectors();
    void generate_hci_queue();
    void generate_orb_lut();
    double get_abs_hamiltonian_by_pqrs(
        const int, const int, const int, const int);
};

}

#endif // HCI_HEG_SOLVER_H_
