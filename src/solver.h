#ifndef HCI_SOLVER_H_
#define HCI_SOLVER_H_

#include "det.h"
#include "wavefunction.h"

namespace hci {

class Solver {
  public:
    void solve();
  protected:
    virtual Wavefunction find_connected_dets(const Det&, const double) = 0;
    virtual double get_hamiltonian_elem(
        const Det&, const Det&,
        const int n_pq_up = -1, const int n_pq_dn = -1,
        const int n_rs_up = -1, const int n_rs_dn = -1) const = 0;
    void pt_det(const double);
    virtual void setup() = 0;

    double max_abs_H;
    int max_n_rs_pairs;
    int n_elecs;
    int n_orbs;
    int n_up;
    int n_dn;
    hci::Wavefunction wf;
    double var_energy;
    double pt_energy;
};

}

#endif // HCI_SOLVER_H_
