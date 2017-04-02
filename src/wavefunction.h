#ifndef HCI_WAVEFUNCTION_H_
#define HCI_WAVEFUNCTION_H_

#include <string>
#include "det.h"

namespace hci {

class Wavefunction {
  public:
    Wavefunction();
    void load(const std::string&, const int);
    void append_det(const Det&, const double coef = 0.0);
    const hci::Det& get_det(const int) const;
    hci::Det& get_det(const int);
    void set_det(const int, const Det&);
    double get_coef(const int) const;
    void set_coef(const int, const double);

    int n; // Number of dets.
  private:
    std::vector<double> coefs;
    std::vector<hci::Det> dets;
};

}

#endif // HCI_WAVEFUNCTION_H_
