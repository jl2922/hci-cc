#ifndef HCI_WAVEFUNCTION_H_
#define HCI_WAVEFUNCTION_H_

#include <string>
#include "det.h"

namespace hci {

class Wavefunction {
  public:
    Wavefunction();
    void load(std::string);
    const hci::Det& get_det(int) const;
    double get_coef(int) const;

    int n; // Number of dets.
  private:
    std::vector<double> coefs;
    std::vector<hci::Det> dets;
};

}

#endif // HCI_WAVEFUNCTION_H_
