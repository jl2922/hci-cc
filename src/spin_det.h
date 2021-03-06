#ifndef HCI_SPIN_DET_H_
#define HCI_SPIN_DET_H_

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/serialization/access.hpp>
#include <cstddef>
#include <forward_list>
#include <iostream>
#include <vector>

namespace hci {

typedef uint32_t BitsBlock;
typedef uint32_t EncodeBlock;
typedef std::vector<EncodeBlock> EncodeType;

class SpinDet {
  friend class boost::serialization::access;
  friend std::ostream &operator<<(std::ostream &, const SpinDet &);

 public:
  SpinDet(){};
  SpinDet(const int n_orb);
  SpinDet(const SpinDet &spinDet);
  friend bool operator==(const SpinDet &, const SpinDet &);
  SpinDet &operator=(const SpinDet &);
  int get_n_elecs() const;
  std::vector<int> get_elec_orbs(const int n_elecs = -1) const;
  int get_n_orbs() const;
  void resize(const int);
  SpinDet &set_orb(const int, const bool occ = true);
  bool get_orb(const int) const;
  void from_eor(const SpinDet &, const SpinDet &);
  bool is_zero() const;
  void print();
  EncodeType encode() const;
  void decode(const EncodeType &, const int n_orbs);
  friend std::size_t hash_value(const SpinDet &);

 private:
  boost::dynamic_bitset<BitsBlock> orbs;
};

bool operator==(const SpinDet &, const SpinDet &);
std::size_t hash_value(const SpinDet &);

std::ostream &operator<<(std::ostream &, const SpinDet &);
}

#endif  // HCI_SPIN_DET_H_
