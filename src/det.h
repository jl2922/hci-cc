#ifndef HCI_DEG_H_
#define HCI_DEG_H_

#include <cstddef>
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
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
    std::vector<BitsBlock> as_vector() const;
    void from_vector(const std::vector<BitsBlock>&, const int);
    hci::SpinDet up;
    hci::SpinDet dn;
  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      ar & up;
      ar & dn;
    }
};

bool operator==(const Det&, const Det&);
std::size_t hash_value(const Det&);

}

#endif // HCI_DEG_H_
