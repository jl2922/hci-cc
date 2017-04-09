#include "det.h"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/functional/hash.hpp>
#include <cstddef>
#include <iostream>

#ifdef BIG_SPIN_DET
#include "big_spin_det.h"
#define SpinDet BigSpinDet
#else
#include "spin_det.h"
#endif

namespace hci {

Det::Det(const int n_orbs) { resize(n_orbs); }

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

bool Det::is_zero() const { return up.is_zero() && dn.is_zero(); }

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

std::pair<EncodeType, EncodeType> Det::encode() const {
  std::pair<EncodeType, EncodeType> pair;
  pair.first = up.encode();
  pair.second = dn.encode();
  return pair;
}

void Det::decode(const std::pair<EncodeType, EncodeType>& pair, const int n_orbs) {
  up.decode(pair.first, n_orbs);
  dn.decode(pair.second, n_orbs);
}

std::size_t hash_value(const Det& det) {
  std::size_t seed = hash_value(det.up);
  boost::hash_combine(seed, hash_value(det.dn));
  return seed;
}
}
