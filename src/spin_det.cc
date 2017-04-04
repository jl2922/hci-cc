#include "spin_det.h"

#include <cstddef>
#include <forward_list>
#include <iostream>
#include <vector>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

namespace hci {

SpinDet::SpinDet(const int n_orbs) {
  resize(n_orbs); // Fill with zeroes.
}

SpinDet::SpinDet(const SpinDet& spin_det) {
  orbs = spin_det.orbs;
}

bool operator==(const SpinDet& lhs, const SpinDet& rhs) {
  return lhs.orbs == rhs.orbs;
}

SpinDet& SpinDet::operator=(const SpinDet& rhs) {
  orbs = rhs.orbs;
  return *this;
}

void SpinDet::from_eor(const SpinDet& lhs, const SpinDet& rhs) {
  orbs = lhs.orbs ^ rhs.orbs;
}

int SpinDet::get_n_elecs() const {
  return orbs.count();
}

std::vector<int> SpinDet::get_elec_orbs(int n_elecs) const {
  std::vector<int> elec_orbs;
  if (n_elecs > 0) elec_orbs.reserve(n_elecs);
  std::size_t pos = orbs.find_first();
  while(pos != boost::dynamic_bitset<>::npos) {
    elec_orbs.push_back(static_cast<int>(pos));
    pos = orbs.find_next(pos);
  }
  return elec_orbs;
}

int SpinDet::get_n_orbs() const {
  return orbs.size();
}

bool SpinDet::is_zero() const {
  return orbs.none();
}

void SpinDet::resize(const int n_orbs) {
  orbs.resize(n_orbs);
}

SpinDet& SpinDet::set_orb(const int orb, const bool occ) {
  orbs.set(orb, occ);
  return *this;
}

bool SpinDet::get_orb(const int orb) const {
  return orbs.test(orb);
}

std::size_t hash_value(const SpinDet& spin_det) {
  std::size_t seed = 0;
  std::size_t pos = spin_det.orbs.find_first(); 
  while (pos != boost::dynamic_bitset<>::npos) {
    boost::hash_combine(seed, pos);
    pos = spin_det.orbs.find_next(pos);
  }
  return seed;
}

}
