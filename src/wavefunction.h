#ifndef HCI_WAVEFUNCTION_H_
#define HCI_WAVEFUNCTION_H_

#include <list>
#include <string>
#include "det.h"

namespace hci {

class Wavefunction {
  public:
    Wavefunction();
    void load(const std::string&, const int);
    Det& append_det(const Det&, const double coef = 0.0);
    const std::list<Det>& get_dets() const;
    const std::list<double>& get_coefs() const;
    void freeze();

    int n; // Number of dets.
  private:
    std::list<double> coefs;
    std::list<hci::Det> dets;
};

}

#endif // HCI_WAVEFUNCTION_H_
