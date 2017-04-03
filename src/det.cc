#include "det.h"

#include <cstddef>
#include <iostream>
#include <boost/functional/hash.hpp>
#include "spin_det.h"

namespace hci {

Det::Det(const int n_orbs) {
  resize(n_orbs);
}
  
Det::Det(const Det& det) {
  up = det.up;
  dn = det.dn;
}

bool operator==(const Det& lhs, const Det& rhs) {
  return lhs.up == rhs.up && lhs.dn == rhs.dn;
}

Det& Det::operator=(const Det& det) {
  up = det.up;
  dn = det.dn;
  return *this;
}

void Det::from_eor(const Det& lhs, const Det& rhs) {
  up.from_eor(lhs.up, rhs.up);
  dn.from_eor(lhs.dn, rhs.dn);
}

bool Det::get_orb(const int orb, const int n_orbs) const {
  if (orb < n_orbs) {
    return up.get_orb(orb);
  } else {
    return dn.get_orb(orb - n_orbs);
  }
}

Det& Det::set_orb(const int orb, const int n_orbs, const bool occ) {
  if (orb < n_orbs) {
    up.set_orb(orb, occ);
  } else {
    dn.set_orb(orb - n_orbs, occ);
  }
  return *this;
}

void Det::resize(const int n_orbs) {
  up.resize(n_orbs);
  dn.resize(n_orbs);
}

std::size_t hash_value(const Det& det) {
  std::size_t seed = hash_value(det.up);
  boost::hash_combine(seed, hash_value(det.dn));
  return seed;
}



}
