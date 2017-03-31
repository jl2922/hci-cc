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
    int get_n_elec() const;
    int* get_elec_orbs(int n_elec = -1) const;
    void resize(const int);
    void set_orb(const int, const bool);
    bool get_orb(const int) const;
    void from_eor(const SpinDet&, const SpinDet&);
    void print();
    friend std::size_t hash_value(const SpinDet&);
  private:
    boost::dynamic_bitset<> orbs;
};

bool operator==(const SpinDet&, const SpinDet&);
std::size_t hash_value(const SpinDet&);

}

#endif // HCI_SPIN_DET_H_
