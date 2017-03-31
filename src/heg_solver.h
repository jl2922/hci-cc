#ifndef HCI_HEG_SOLVER_H_
#define HCI_HEG_SOLVER_H_

#include "det.h"
#include "solver.h"
#include "wavefunction.h"

namespace hci {

struct HEGData {
  double cell_length;
  double k_unit;
  std::vector<std::array<int, 3>> k_vectors;
  int n_diff;
  int n_diff_offset;
  int n_max;
  double r_cutoff;
  double r_s;
};

class HEGSolver: public Solver {
  public:
  protected:
    Wavefunction find_connected_dets(const Det&, const double);
    double get_hamiltonian_elem(const Det&, const Det&);
    void setup();
    HEGData heg;
  private:
    int get_gamma_exp(const SpinDet&, const int*, const int*, const int);
    int generate_k_vectors();
};

}

#endif // HCI_HEG_SOLVER_H_
