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
    void pt_det(const double);
    virtual void setup() = 0;
    hci::Wavefunction wf;
};

}

#endif // HCI_SOLVER_H_