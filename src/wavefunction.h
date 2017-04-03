#ifndef HCI_WAVEFUNCTION_H_
#define HCI_WAVEFUNCTION_H_

#include <list>
#include "det.h"

namespace hci {

class Wavefunction {
  public:
    Wavefunction();
    Det& append_det(const Det&, const double coef = 0.0);
    const std::list<Det>& get_dets() const;
    const std::list<double>& get_coefs() const;
    int size() const;
  private:
    std::list<double> coefs;
    std::list<hci::Det> dets;
    int n; // Number of dets.
};

}

#endif // HCI_WAVEFUNCTION_H_
