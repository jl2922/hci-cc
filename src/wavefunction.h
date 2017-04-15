#ifndef HCI_WAVEFUNCTION_H_
#define HCI_WAVEFUNCTION_H_

#include <list>
#include "det.h"

namespace hci {

class Wavefunction {
 public:
  Wavefunction() { n = 0; }
  Det& append_det(const Det& det, const double coef = 0.0) {
    n++;
    dets.push_back(det);
    coefs.push_back(coef);
    return dets.back();
  }
  const std::list<Det>& get_dets() const { return dets; }
  const std::list<double>& get_coefs() const { return coefs; }
  void clear() {
    dets.clear();
    coefs.clear();
    n = 0;
  }
  int size() const { return n; }

 private:
  std::list<double> coefs;
  std::list<Det> dets;
  int n;  // Number of dets.
};
}

#endif  // HCI_WAVEFUNCTION_H_
