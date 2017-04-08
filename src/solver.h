#ifndef HCI_SOLVER_H_
#define HCI_SOLVER_H_

#include <boost/mpi.hpp>
#include <list>
#include <string>
#include "det.h"
#include "wavefunction.h"

namespace hci {

struct MpiData {
  int id;
  int n;
  boost::mpi::environment env;
  boost::mpi::communicator world;
};

class Solver {
 public:
  Solver();
  void solve();

 protected:
  virtual std::list<Det> find_connected_dets(const Det &, const double) = 0;
  virtual double get_hamiltonian_elem(
      const Det &,
      const Det &,
      const int n_up = -1,
      const int n_dn = -1) const = 0;
  void load_wavefunction(const std::string &);
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
  MpiData mpi;
};
}

#endif  // HCI_SOLVER_H_
