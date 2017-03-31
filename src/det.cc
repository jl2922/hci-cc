#include "det.h"

#include <cstddef>
#include <boost/functional/hash.hpp>
#include "spin_det.h"

namespace hci {
  
Det::Det(const Det& det) {
  this->up = det.up;
  this->dn = det.dn;
}

bool operator==(const Det& det1, const Det& det2) {
  return det1.up == det2.up && det1.dn == det2.dn;
}

Det& Det::operator=(const Det& det) {
  this->up = det.up;
  this->dn = det.dn;
  return *this;
}

void Det::from_eor(const Det& det_a, const Det& det_b) {
  this->up.from_eor(det_a.up, det_b.up);
  this->dn.from_eor(det_a.dn, det_b.dn);
}

void Det::resize(const int n_orb) {
  this->up.resize(n_orb);
  this->dn.resize(n_orb);
}

std::size_t hash_value(const Det& det) {
  std::size_t seed = hash_value(det.up);
  boost::hash_combine(seed, hash_value(det.dn));
  return seed;
}

}
