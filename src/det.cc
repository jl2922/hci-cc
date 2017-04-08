#include "det.h"

#include <cstddef>
#include <iostream>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
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

bool Det::is_zero() const {
  return up.is_zero() && dn.is_zero();
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

std::vector<BitsBlock> Det::as_vector() const {
  std::vector<BitsBlock> vec;
  vec.reserve(up.orbs.num_blocks() + dn.orbs.num_blocks());
  to_block_range(up.orbs, std::back_insert_iterator<std::vector<BitsBlock>>(vec));
  to_block_range(dn.orbs, std::back_insert_iterator<std::vector<BitsBlock>>(vec));
  // for (const auto& x: vec) std::cout << x << " ";
  // std::cout << std::endl;
  // Det tmp;
  // tmp.from_vector(vec, 81);

  return vec;
}

void Det::from_vector(const std::vector<BitsBlock>& vec, const int n_orbs) {
  const int num_blocks = vec.size();

  up.orbs.clear();
  up.orbs.reserve(n_orbs);
  for (int i = 0; i < num_blocks / 2; i++) up.orbs.append(vec[i]);
  up.orbs.resize(n_orbs);
  // std::cout << "up:" << up.orbs << std::endl;

  dn.orbs.clear();
  dn.orbs.reserve(n_orbs);
  for (int i = num_blocks / 2; i < num_blocks; i++) dn.orbs.append(vec[i]);
  dn.orbs.resize(n_orbs);
  // std::cout << "dn:" << dn.orbs << std::endl;

}

std::size_t hash_value(const Det& det) {
  std::size_t seed = hash_value(det.up);
  boost::hash_combine(seed, hash_value(det.dn));
  return seed;
}



}
