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
    virtual double get_hamiltonian_elem(const Det&, const Det&) = 0;
    void pt_det(const double);
    virtual void setup() = 0;
    double max_abs_H;
    int n_dn;
    int n_elec;
    int n_orb;
    int n_up;
    hci::Wavefunction wf;
};

}

#endif // HCI_SOLVER_H_
