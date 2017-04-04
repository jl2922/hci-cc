#include <iostream>

#include "parallel.h"
#include "heg_solver.h"

int main(int argc, char **argv) {
  hci::Parallel::getInstance().init(argc, argv);
  hci::HEGSolver heg_solver;
  heg_solver.solve();
  hci::Parallel::getInstance().finish();
}
