#ifndef HCI_SPIN_DET_H_
#define HCI_SPIN_DET_H_

#include <cstddef>
#include <iostream>
#include <vector>
#include <boost/dynamic_bitset.hpp>

namespace hci {

class SpinDet {
  public:
    SpinDet() { };
    SpinDet(const int n_orb);
    SpinDet(const SpinDet& spinDet);
    friend bool operator==(const SpinDet&, const SpinDet&);
    SpinDet& operator=(const SpinDet&);
    int get_n_elecs() const;
    std::vector<int> get_elec_orbs(const int n_elecs = -1) const;
    int get_n_orbs() const;
    void resize(const int);
    SpinDet& set_orb(const int, const bool occ = true);
    bool get_orb(const int) const;
    void from_eor(const SpinDet&, const SpinDet&);
    void print();
    friend std::size_t hash_value(const SpinDet&);
//    friend std::ostream& operator<<(std::ostream&, const SpinDet&);
  private:
    boost::dynamic_bitset<> orbs;
};

bool operator==(const SpinDet&, const SpinDet&);
std::size_t hash_value(const SpinDet&);
//std::ostream& operator<<(std::ostream&, const SpinDet&);

}

#endif // HCI_SPIN_DET_H_
