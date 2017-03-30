#ifndef HCI_DEG_H_
#define HCI_DEG_H_

#include <cstddef>
#include "spin_det.h"

namespace hci {

class Det {
  public:
    Det() { };
    Det(const Det&);
    Det& operator=(const Det&);
    void resize(const int);
    hci::SpinDet up;
    hci::SpinDet dn;
  private:
};

bool operator==(const Det&, const Det&);
std::size_t hash_value(const Det&);

}

#endif // HCI_DEG_H_