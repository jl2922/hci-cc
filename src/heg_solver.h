#ifndef HCI_HEG_SOLVER_H_
#define HCI_HEG_SOLVER_H_

#include "det.h"
#include "solver.h"
#include "wavefunction.h"

namespace hci {

class HEGSolver: public Solver {
  public:
  protected:
    Wavefunction find_connected_dets(const Det& det, const double eps);
    void setup();
};

}

#endif // HCI_HEG_SOLVER_H_