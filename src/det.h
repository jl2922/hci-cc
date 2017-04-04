#ifndef HCI_DEG_H_
#define HCI_DEG_H_

#include <cstddef>
#include <iostream>
#include "spin_det.h"

namespace hci {

class Det {
  public:
    Det() { };
    Det(const int);
    Det(const Det&);
    Det& operator=(const Det&);
    bool get_orb(const int, const int) const;
    void from_eor(const Det&, const Det&);
    bool is_zero() const;
    void resize(const int);
    Det& set_orb(const int, const int, const bool occ = true);
    
    hci::SpinDet up;
    hci::SpinDet dn;
  private:
};

bool operator==(const Det&, const Det&);
std::size_t hash_value(const Det&);

}

#endif // HCI_DEG_H_
