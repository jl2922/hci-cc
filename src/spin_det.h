#ifndef HCI_SPIN_DET_H_
#define HCI_SPIN_DET_H_

#include <cstddef>
#include <boost/dynamic_bitset.hpp>

namespace hci {

class SpinDet {
  public:
    SpinDet() { };
    SpinDet(const int n_orb);
    SpinDet(const SpinDet& spinDet);
    friend bool operator==(const SpinDet&, const SpinDet&);
    SpinDet& operator=(const SpinDet&);
    void resize(const int n_orb);
    void set_orb(const int orb, const bool occ = true);
    void print();
    friend std::size_t hash_value(const SpinDet&);
  private:
    boost::dynamic_bitset<> orbs;
};

bool operator==(const SpinDet&, const SpinDet&);
std::size_t hash_value(const SpinDet&);

}

#endif // HCI_SPIN_DET_H_