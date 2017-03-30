#include "spin_det.h"

#include <cstddef>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

namespace hci {

SpinDet::SpinDet(const int n_orb) {
  this->resize(n_orb);
}

SpinDet::SpinDet(const SpinDet& spin_det) {
  this->orbs = spin_det.orbs;
}

bool operator==(const SpinDet& spin_det1, const SpinDet& spin_det2) {
  return spin_det1.orbs == spin_det2.orbs;
}

SpinDet& SpinDet::operator=(const SpinDet& spin_det) {
  this->orbs = spin_det.orbs;
  return *this;
}

void SpinDet::resize(const int n_orb) {
  this->orbs.resize(n_orb);
}

void SpinDet::set_orb(const int orb, const bool occ) {
  this->orbs.set(orb, occ);
}

void SpinDet::print() {
  std::cout << this->orbs << std::endl;
}

std::size_t hash_value(const SpinDet& spin_det) {
  std::size_t seed = 0;
  std::size_t pos = spin_det.orbs.find_first(); 
  while (pos != boost::dynamic_bitset<>::npos) {
    boost::hash_combine(seed, spin_det.orbs.find_next(pos));
    pos = spin_det.orbs.find_next(pos);
  }
  return seed;
}

}
