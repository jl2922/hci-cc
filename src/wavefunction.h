#ifndef HCI_WAVEFUNCTION_H_
#define HCI_WAVEFUNCTION_H_

#include <string>
#include "det.h"

namespace hci {

class Wavefunction {
  public:
    Wavefunction();
    void load(std::string);
    hci::Det& get_det(int);
    double get_coef(int);
    int n;
  private:
    double* coefs;
    hci::Det* dets;
};

}

#endif // HCI_WAVEFUNCTION_H_